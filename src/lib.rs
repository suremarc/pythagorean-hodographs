#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

use alloc::vec::Vec;
use core::ops::{Add, Mul};
use glam::f32::*;
use itertools::Itertools;
#[cfg(not(feature = "std"))]
#[allow(unused_imports)]
use num_traits::Float;
use ordered_float::OrderedFloat;

/// A twice-differentiable parametric curve in 3-space.
pub trait Curve3 {
    /// The value of this curve.
    fn p(&self, u: f32) -> Vec3;
    /// The derivative of this curve.
    fn dp(&self, u: f32) -> Vec3;
    /// The second derivative of this curve.
    fn d2p(&self, u: f32) -> Vec3;

    /// The speed of this curve.
    /// Pythagorean hodographs have speed expressible as a simple polynomial, which makes
    /// certain quantities (such as arc length) efficient to compute.
    fn speed(&self, u: f32) -> f32;

    // The arc length of this curve between 0 and u.
    fn arc_length(&self, u: f32) -> f32;

    fn tangent(&self, u: f32) -> Vec3 {
        self.dp(u) / self.speed(u)
    }
}

/// A parametric moving frame in 3-space.
pub trait Frame: Curve3 {
    fn quat_unnormalized(&self, u: f32) -> Quat;

    // A quaternion representation of this moving frame's rotational component.
    fn quat(&self, u: f32) -> Quat {
        self.quat_unnormalized(u).normalize()
    }

    /// The value of moving frame in 3-space, containing both translation and rotation components.
    fn frame(&self, u: f32) -> Affine3A {
        let a = self.quat_unnormalized(u);

        let ai = a.mul_vec3(vec3(1., 0., 0.)) / a.length_squared();
        let aj = a.mul_vec3(vec3(0., 1., 0.)) / a.length_squared();
        let ak = a.mul_vec3(vec3(0., 0., 1.)) / a.length_squared();

        Affine3A::from_mat3_translation(Mat3::from_cols(ai, aj, ak), self.p(u))
    }
}

/// A continuous curve built up from piecewise polynomials.
#[derive(Debug, Clone)]
pub struct Spline<T: Curve3> {
    pub segments: Vec<T>,
}

impl<T: Curve3> Spline<T> {
    fn normalize(&self, mut u: f32) -> (f32, usize) {
        u *= self.segments.len() as f32;
        let i = u.floor().clamp(0., self.segments.len() as f32 - 1.);

        (u - i, i as usize)
    }
}

impl<T: Curve3> Curve3 for Spline<T> {
    #[inline]
    fn p(&self, u: f32) -> Vec3 {
        let (u, i) = self.normalize(u);
        self.segments[i].p(u)
    }

    #[inline]
    fn dp(&self, u: f32) -> Vec3 {
        let (u, i) = self.normalize(u);
        self.segments[i].dp(u) * self.segments.len() as f32
    }

    #[inline]
    fn d2p(&self, u: f32) -> Vec3 {
        let (u, i) = self.normalize(u);
        self.segments[i].d2p(u) * self.segments.len().pow(2) as f32
    }

    #[inline]
    fn speed(&self, u: f32) -> f32 {
        let (u, i) = self.normalize(u);
        self.segments[i].speed(u) * self.segments.len() as f32
    }

    // TODO: make this more efficient
    fn arc_length(&self, u: f32) -> f32 {
        let (u, i) = self.normalize(u);
        self.segments[0..i]
            .iter()
            .map(|h| h.arc_length(1.0))
            .sum::<f32>()
            + self.segments[i].arc_length(u)
    }
}

impl<T: Curve3> Frame for Spline<T>
where
    T: Frame,
{
    #[inline]
    fn quat_unnormalized(&self, u: f32) -> Quat {
        let (u, i) = self.normalize(u);
        self.segments[i].quat_unnormalized(u)
    }
}

#[derive(Debug, Clone)]
pub struct QuinticPHCurve {
    data: QuinticPHData,
    hermite: HermiteQuintic,
}

impl QuinticPHCurve {
    /// Construct a new quintic PH curve given C1 data.
    /// pi and pf are the initial and final positions, while
    /// di and df are the initial and final derivatives.
    fn new(pi: Vec3, pf: Vec3, di: Vec3, df: Vec3) -> Self {
        let dp = pf - pi;
        let (q, _) = inverse_solve_quat(di + df);
        let q_normalized = q.normalize();
        let q_normalized_inv = q_normalized.conjugate();
        let (tdi, tdf) = (q_normalized_inv.mul_vec3(di), q_normalized_inv.mul_vec3(df));

        let (ta0, _) = inverse_solve_quat(tdi);
        let (ta2, _) = inverse_solve_quat(tdf);

        [(ta0, ta2), (-ta0, ta2), (ta0, -ta2), (-ta0, -ta2)]
            .into_iter()
            .map(|(ta0, ta2)| {
                const F: Quat = Quat::from_array([0., 0., -1., 0.]); // forward
                let c = 120. * q_normalized_inv.mul_vec3(dp) - 15. * (tdi + tdf)
                    + 5. * (ta0 * F * ta2.conjugate() + ta2 * F * ta0.conjugate()).xyz();
                let (r, _) = inverse_solve_quat(c);
                let ta1 = (ta0 + ta2) * -0.75 + r * 0.25;

                (ta0, ta1, ta2)
            })
            .map(|(ta0, ta1, ta2)| (q_normalized * ta0, q_normalized * ta1, q_normalized * ta2))
            .map(|(a0, a1, a2)| QuinticPHData { a0, a1, a2, pi })
            .map(|data| Self {
                data,
                hermite: data.curve(),
            })
            .min_by_key(|h| OrderedFloat(h.hermite.elastic_bending_energy()))
            .unwrap()
    }
}

impl Curve3 for QuinticPHCurve {
    #[inline]
    fn p(&self, u: f32) -> Vec3 {
        self.hermite.p(u)
    }

    #[inline]
    fn dp(&self, u: f32) -> Vec3 {
        self.hermite.dp(u)
    }

    #[inline]
    fn d2p(&self, u: f32) -> Vec3 {
        self.hermite.d2p(u)
    }

    #[inline]
    fn speed(&self, u: f32) -> f32 {
        self.hermite.speed(u)
    }

    #[inline]
    fn arc_length(&self, u: f32) -> f32 {
        self.hermite.arc_length(u)
    }
}

impl Frame for QuinticPHCurve {
    /// The quaternion associated with the Euler-Rodrigues frame for this PH curve.
    /// Note that a spline consisting of Euler-Rodrigues frames is not necessarily continuous;
    /// in the spline implementation we adjust the frames such that they are continuous.
    fn quat_unnormalized(&self, u: f32) -> Quat {
        self.data.a0 * (1. - u).powi(2)
            + self.data.a1 * u * (1. - u) * 2.
            + self.data.a2 * u.powi(2)
    }
}

impl Spline<QuinticPHCurve> {
    /// Construct a quintic PH spline from positional data, using the Catmull-Rom method
    /// for the differential data.
    pub fn catmull_rom(pts: &[Vec3]) -> Self {
        Self {
            segments: pts
                .iter()
                .copied()
                .tuple_windows()
                .map(move |(p0, p1, p2)| (p1, 0.5 * (p2 - p0)))
                .tuple_windows()
                .map(|((p0, d0), (p1, d1))| QuinticPHCurve::new(p0, p1, d0, d1))
                .collect(),
        }
        .make_frames_continuous()
    }

    /// Construct a quintic PH spline given C1 interpolation data (position and derivative pairs).
    pub fn hermite(data: &[(Vec3, Vec3)]) -> Self {
        Self {
            segments: data
                .iter()
                .copied()
                .tuple_windows()
                .map(|((p0, d0), (p1, d1))| QuinticPHCurve::new(p0, p1, d0, d1))
                .collect(),
        }
        .make_frames_continuous()
    }

    fn make_frames_continuous(mut self) -> Self {
        for i in 1..self.segments.len() {
            let delta_a = self.segments[i].data.a0.conjugate() * self.segments[i - 1].data.a2;
            let data = &mut self.segments[i].data;
            for a in [&mut data.a0, &mut data.a1, &mut data.a2] {
                *a *= delta_a;
            }
        }

        self
    }
}

#[derive(Debug, Clone, Copy)]
struct QuinticPHData {
    a0: Quat,
    a1: Quat,
    a2: Quat,
    pi: Vec3,
}

impl QuinticPHData {
    fn curve(&self) -> HermiteQuintic {
        let a0len2 = self.a0.length_squared();
        let a1len2 = self.a1.length_squared();
        let a2len2 = self.a2.length_squared();
        let a0a1 = 2. * self.a0.dot(self.a1);
        let a0a2 = 2. * self.a0.dot(self.a2);
        let a1a2 = 2. * self.a1.dot(self.a2);

        let w0 = a0len2;
        let w1 = 0.5 * a0a1;
        let w2 = (1. / 6.) * (a0a2 + 4. * a1len2);
        let w3 = 0.5 * a1a2;
        let w4 = a2len2;

        const F: Quat = Quat::from_array([0., 0., -1., 0.]);

        let a0fa0 = (self.a0 * F * self.a0.conjugate()).xyz();
        let a1fa1 = (self.a1 * F * self.a1.conjugate()).xyz();
        let a2fa2 = (self.a2 * F * self.a2.conjugate()).xyz();
        let a0fa1 = (self.a0 * F * self.a1.conjugate() + self.a1 * F * self.a0.conjugate()).xyz();
        let a0fa2 = (self.a0 * F * self.a2.conjugate() + self.a2 * F * self.a0.conjugate()).xyz();
        let a1fa2 = (self.a1 * F * self.a2.conjugate() + self.a2 * F * self.a1.conjugate()).xyz();

        let wt0 = a0fa0;
        let wt1 = 0.5 * a0fa1;
        let wt2 = (1. / 6.) * (a0fa2 + 4. * a1fa1);
        let wt3 = 0.5 * a1fa2;
        let wt4 = a2fa2;

        HermiteQuintic {
            pi: self.pi,
            weights: [w0, w1, w2, w3, w4],
            weighted_tangents: [wt0, wt1, wt2, wt3, wt4],
        }
    }
}

/// Hermite quintic polynomial obtained from the construction of a quintic PH spline.
#[derive(Debug, Clone)]
struct HermiteQuintic {
    pi: Vec3,
    weights: [f32; 5],
    weighted_tangents: [Vec3; 5],
}

impl HermiteQuintic {
    // TODO: investigate tradeoffs in precision and speed (number of subdivisions)
    // exact methods are available but are so complex to compute that
    // it's equally as efficient to do a Riemann sum with 500 subdivisions
    fn elastic_bending_energy(&self) -> f32 {
        (0..25)
            .map(|x| x as f32 * 0.04)
            .map(|u| self.dp(u).cross(self.d2p(u)).length_squared() / self.speed(u).powi(5) * 0.04)
            .sum()
    }
}

impl Curve3 for HermiteQuintic {
    fn p(&self, u: f32) -> Vec3 {
        self.pi + hermite_quintic_polynomial_integral(self.weighted_tangents, u)
            - hermite_quintic_polynomial_integral(self.weighted_tangents, 0.)
    }

    fn dp(&self, u: f32) -> Vec3 {
        hermite_quintic_polynomial(self.weighted_tangents, u)
    }

    fn d2p(&self, u: f32) -> Vec3 {
        hermite_quintic_polynomial_derivative(self.weighted_tangents, u)
    }

    fn speed(&self, u: f32) -> f32 {
        hermite_quintic_polynomial(self.weights, u)
    }

    fn arc_length(&self, u: f32) -> f32 {
        hermite_quintic_polynomial_integral(self.weights, u)
    }
}

/// Generates a basis for the one-parameter family of solutions to -AkA* = c.
/// This function returns the quaternions A1, A2 such that any solution
/// to -AkA* = c can be generated by the solution A1 cos t + A2 sin t.
// TODO: use an approximation near nu = 1.0
fn inverse_solve_quat(c: Vec3) -> (Quat, Quat) {
    let [lambda, mu, nu] = c.normalize().to_array();
    let r = (0.5 * (1. - nu) * c.length()).sqrt();
    (
        quat(-lambda / (1. - nu), -mu / (1. - nu), 1., 0.) * r,
        quat(mu / (1. - nu), -lambda / (1. - nu), 0., 1.) * r,
    )
}

fn hermite_quintic_polynomial<T: Copy + Add<Output = T> + Mul<f32, Output = T>>(
    coefs: [T; 5],
    t: f32,
) -> T {
    coefs[0] * (1. - t).powi(4)
        + coefs[1] * (1. - t).powi(3) * t * 4.
        + coefs[2] * ((1. - t) * t).powi(2) * 6.
        + coefs[3] * (1. - t) * t.powi(3) * 4.
        + coefs[4] * t.powi(4)
}

fn hermite_quintic_polynomial_integral<T: Copy + Add<Output = T> + Mul<f32, Output = T>>(
    coefs: [T; 5],
    t: f32,
) -> T {
    coefs[0] * 0.2 * (t - 1.).powi(5)
        + coefs[1] * (-0.2 * t.powi(5) + 0.75 * t.powi(4) - t.powi(3) + 0.5 * t.powi(2)) * 4.
        + coefs[2] * (0.2 * t.powi(5) - 0.5 * t.powi(4) + t.powi(3) / 3.) * 6.
        + coefs[3] * (0.25 * t.powi(4) - 0.2 * t.powi(5)) * 4.
        + coefs[4] * 0.2 * t.powi(5)
}

fn hermite_quintic_polynomial_derivative<T: Copy + Add<Output = T> + Mul<f32, Output = T>>(
    coefs: [T; 5],
    t: f32,
) -> T {
    coefs[0] * -4. * (1. - t).powi(3)
        + coefs[1] * (1. - 4. * t) * (1. - t).powi(2) * 4.
        + coefs[2] * 2. * (t - 1.) * t * (2. * t - 1.) * 6.
        + coefs[3] * (3. - 4. * t) * t.powi(2) * 4.
        + coefs[4] * 4. * t.powi(3)
}

#[cfg(test)]
mod tests {}
