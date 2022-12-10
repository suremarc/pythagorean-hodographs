use std::ops::{Add, Mul};

use glam::f32::*;

use itertools::Itertools;

// A smoothly varying curve in 3-space.
pub trait Curve {
    fn p(&self, u: f32) -> Vec3;
    fn speed(&self, u: f32) -> f32;
}

// A single-parameter moving frame in 3-space.
pub trait Frame {
    fn frame(&self, u: f32) -> Affine3A;
}

// A continuous curve built up from piecewise polynomials.
pub struct Spline<T: Curve> {
    segments: Vec<T>,
}

impl<T: Curve> Spline<T> {
    #[inline]
    fn normalize(&self, mut u: f32) -> (f32, usize) {
        u *= self.segments.len() as f32;
        let i = u.floor().clamp(0., self.segments.len() as f32 - 1.);

        (u - i, i as usize)
    }
}

impl<T: Curve> Curve for Spline<T> {
    #[inline]
    fn p(&self, u: f32) -> Vec3 {
        let (u, i) = self.normalize(u);
        self.segments[i].p(u)
    }

    #[inline]
    fn speed(&self, u: f32) -> f32 {
        let (u, i) = self.normalize(u);
        self.segments[i].speed(u) * self.segments.len() as f32
    }
}

impl<T: Curve> Frame for Spline<T>
where
    T: Frame,
{
    fn frame(&self, u: f32) -> Affine3A {
        let (u, i) = self.normalize(u);
        self.segments[i].frame(u)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct QuinticPHCurve {
    a0: Quat,
    a1: Quat,
    a2: Quat,
    pi: Vec3,
}

impl QuinticPHCurve {
    pub fn new(pi: Vec3, pf: Vec3, di: Vec3, df: Vec3) -> QuinticPHCurve {
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
                const I: Quat = Quat::from_array([1., 0., 0., 0.]);
                let c = 120. * q_normalized_inv.mul_vec3(dp) - 15. * (tdi + tdf)
                    + 5. * (ta0 * I * ta2.conjugate() + ta2 * I * ta0.conjugate()).xyz();
                let [l, mu, nu] = c.normalize().to_array();
                let ta1 = (ta0 + ta2) * -0.75
                    + Quat::from_xyzw(1., mu / (1. + l), nu / (1. + l), 0.)
                        * 0.25
                        * (0.5 * (1. + l) * c.length()).sqrt();

                (ta0, ta1, ta2)
            })
            .map(|(ta0, ta1, ta2)| (q_normalized * ta0, q_normalized * ta1, q_normalized * ta2))
            .map(|(a0, a1, a2)| QuinticPHCurve { a0, a1, a2, pi })
            .min_by_key(|h| ordered_float::OrderedFloat(h.curve_private().elastic_bending_energy()))
            .unwrap()
    }

    #[inline]
    pub fn euler_rodrigues_frame(&self) -> EulerRodriguesFrame {
        EulerRodriguesFrame {
            data: *self,
            curve: self.curve_private(),
        }
    }

    fn curve_private(&self) -> HermiteQuintic {
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

        const I: Quat = Quat::from_array([1., 0., 0., 0.]);

        let a0ia0 = (self.a0 * I * self.a0.conjugate()).xyz();
        let a1ia1 = (self.a1 * I * self.a1.conjugate()).xyz();
        let a2ia2 = (self.a2 * I * self.a2.conjugate()).xyz();
        let a0ia1 = (self.a0 * I * self.a1.conjugate() + self.a1 * I * self.a0.conjugate()).xyz();
        let a0ia2 = (self.a0 * I * self.a2.conjugate() + self.a2 * I * self.a0.conjugate()).xyz();
        let a1ia2 = (self.a1 * I * self.a2.conjugate() + self.a2 * I * self.a1.conjugate()).xyz();

        let wt0 = a0ia0;
        let wt1 = 0.5 * a0ia1;
        let wt2 = (1. / 6.) * (a0ia2 + 4. * a1ia1);
        let wt3 = 0.5 * a1ia2;
        let wt4 = a2ia2;

        HermiteQuintic {
            pi: self.pi,
            weights: [w0, w1, w2, w3, w4],
            weighted_tangents: [wt0, wt1, wt2, wt3, wt4],
        }
    }

    #[inline]
    pub fn curve(&self) -> impl Curve {
        self.curve_private()
    }
}

pub struct EulerRodriguesFrame {
    data: QuinticPHCurve,
    curve: HermiteQuintic,
}

impl Curve for EulerRodriguesFrame {
    #[inline]
    fn p(&self, u: f32) -> Vec3 {
        self.curve.p(u)
    }

    #[inline]
    fn speed(&self, u: f32) -> f32 {
        self.curve.speed(u)
    }
}

impl Frame for EulerRodriguesFrame {
    fn frame(&self, u: f32) -> Affine3A {
        let a = self.data.a0 * (1. - u).powi(2)
            + self.data.a1 * u * (1. - u) * 2.
            + self.data.a2 * u.powi(2);

        let ai = (a * Quat::from_xyzw(1., 0., 0., 0.) * a.conjugate()).xyz() / a.length_squared();
        let aj = (a * Quat::from_xyzw(0., 1., 0., 0.) * a.conjugate()).xyz() / a.length_squared();
        let ak = (a * Quat::from_xyzw(0., 0., 1., 0.) * a.conjugate()).xyz() / a.length_squared();

        Affine3A::from_mat3_translation(Mat3::from_cols(ai, aj, ak), self.curve.p(u))
    }
}

impl Spline<EulerRodriguesFrame> {
    pub fn catmull_rom(pts: Vec<Vec3>) -> Self {
        let n = pts.len();
        Self {
            segments: pts
                .into_iter()
                .tuple_windows()
                .map(|(p0, p1, p2)| (p1, 0.125 * ((n - 2) as f32) * (p2 - p0)))
                .tuple_windows()
                .map(|((p0, d0), (p1, d1))| {
                    QuinticPHCurve::new(p0, p1, d0, d1).euler_rodrigues_frame()
                })
                .collect(),
        }
    }

    pub fn hermite(data: Vec<[Vec3; 2]>) -> Self {
        Self {
            segments: data
                .into_iter()
                .tuple_windows()
                .map(|([p0, d0], [p1, d1])| {
                    QuinticPHCurve::new(p0, p1, d0, d1).euler_rodrigues_frame()
                })
                .collect(),
        }
    }
}

// Hermite quintic polynomial obtained from the construction of a quintic PH spline.
#[derive(Clone, Debug)]
pub struct HermiteQuintic {
    pi: Vec3,
    weights: [f32; 5],
    weighted_tangents: [Vec3; 5],
}

impl HermiteQuintic {
    #[inline]
    pub fn tangent(&self, u: f32) -> Vec3 {
        hermite_quintic_polynomial(self.weighted_tangents, u)
            / hermite_quintic_polynomial(self.weights, u)
    }

    #[inline]
    pub fn speed(&self, u: f32) -> f32 {
        hermite_quintic_polynomial(self.weights, u)
    }

    #[inline]
    pub fn p(&self, u: f32) -> Vec3 {
        self.pi + hermite_quintic_polynomial_integral(self.weighted_tangents, u)
            - hermite_quintic_polynomial_integral(self.weighted_tangents, 0.)
    }

    #[inline]
    pub fn dp(&self, u: f32) -> Vec3 {
        hermite_quintic_polynomial(self.weighted_tangents, u)
    }

    #[inline]
    pub fn d2p(&self, u: f32) -> Vec3 {
        hermite_quintic_polynomial_derivative(self.weighted_tangents, u)
    }

    pub fn elastic_bending_energy(&self) -> f32 {
        (0..1000)
            .map(|x| x as f32 * 0.001)
            .map(|u| self.curvature_squared(u) * self.speed(u) * 0.001)
            .sum()
    }

    #[inline]
    fn curvature_squared(&self, u: f32) -> f32 {
        self.dp(u).cross(self.d2p(u)).length_squared()
            / hermite_quintic_polynomial(self.weights, u).powi(6)
    }
}

impl Curve for HermiteQuintic {
    fn p(&self, u: f32) -> Vec3 {
        self.p(u)
    }

    fn speed(&self, u: f32) -> f32 {
        self.speed(u)
    }
}

// Generates a basis for the one-parameter family of solutions to AiA* = c.
// This function returns the quaternions A1, A2 such that any solution
// to AiA* = c can be generated by the solution A1 cos t + A2 sin t.
fn inverse_solve_quat(c: Vec3) -> (Quat, Quat) {
    let [lambda, mu, nu] = c.normalize().to_array();
    let r = (0.5 * (1. + lambda) * c.length()).sqrt();
    (
        quat(1., mu / (1. + lambda), nu / (1. + lambda), 0.) * r,
        quat(0., nu / (1. + lambda), -mu / (1. + lambda), -1.) * r,
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
mod tests {
    #[test]
    fn it_works() {}
}
