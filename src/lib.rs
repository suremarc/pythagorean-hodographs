#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;

use alloc::vec::Vec;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Shl};
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
    pub fn new(pi: Vec3, pf: Vec3, di: Vec3, df: Vec3) -> Self {
        let data = QuinticPHData::new(pi, pf, di, df);
        let hermite = data.curve();
        Self { data, hermite }
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

#[derive(Debug, Clone)]
struct QuinticPHData {
    a0: Quat,
    a1: Quat,
    a2: Quat,
    pi: Vec3,
}

impl QuinticPHData {
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
            .min_by_key(|h| OrderedFloat(h.curve().elastic_bending_energy()))
            .unwrap()
    }

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

    fn curvature_squared_numerator(&self) -> Polynomial<14> {
        todo!()
    }

    fn elastic_bending_energy(&self) -> f32 {
        let d0 = self.a0 - self.a1 * 2. + self.a2;
        let d1 = (self.a1 - self.a0) * 2.;
        let d2 = self.a0;
        let (r0, r1) = quaternion_quadratic_roots(d0, d1, d2);

        // assume r0 and r1 are complex for now
        // compute partial fraction decomposition using Gauss-Jordan elimination

        let mut coefs = [Polynomial::<21>::default(); 22];

        let rp = [
            Polynomial([1., -2. * r0.w, r0.length_squared()]),
            Polynomial([1., -2. * r1.w, r1.length_squared()]),
        ];
        let rp5 = rp.map(Polynomial::resize::<21>).map(|x| x * x * x * x * x);

        // index of c_ijk in coefs = 2*(5*i + j) - k
        // i in 0..=1
        // j in 1..=5
        // k in 0..=1
        // coefs[0] = b (the const coefficient)

        coefs[0] = rp5[0] * rp5[1];

        for i in 0..=1 {
            // j=1 is a special case, handle it afterwards
            let mut product = Polynomial::<21>::identity();
            for j in (2..=5).rev() {
                for k in 0..=1 {
                    coefs[2 * (5 * i + j) - k] = (product
                        * -2.
                        * (j - 1) as f32
                        * Polynomial([1., -r0.w]).resize()
                        * rp5[1 - i].resize())
                        << k;
                }

                product *= rp[i].resize();
            }

            for k in 0..=1 {
                coefs[2 * (5 * i + 1) - k] = product << k;
            }
        }

        coefs[21] = self.curvature_squared_numerator().resize();

        let mut system = transpose(coefs.map(|p| p.0));
        gauss_jordan_eliminate(&mut system);

        let solution = system.map(|row| row[21]);

        todo!()
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
    // FIXME: replace this with exact formula for greater precision and speed
    fn elastic_bending_energy(&self) -> f32 {
        (0..1000)
            .map(|x| x as f32 * 0.001)
            .map(|u| self.dp(u).cross(self.d2p(u)).length_squared() / self.speed(u).powi(5) * 0.001)
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

fn quaternion_quadratic_roots(a0: Quat, a1: Quat, a2: Quat) -> (Quat, Quat) {
    let (p, q) = (
        a0.conjugate() / a0.length_squared() * a1,
        a0.conjugate() / a0.length_squared() * a2,
    );

    const EPS: f32 = 0.001;

    let p_prime = quat(p.x, p.y, p.z, 0.);
    let q_prime = q - (p - quat(0., 0., 0., p.w / 2.)) * p.w / 2.;

    if p_prime.length_squared() < 1e-7 {
        // p is real
        let q_im = q - quat(0., 0., 0., q.w);
        if q_im.length_squared() < 1e-7 {
            // both are real; use the regular quadratic formula
            // x^2 + px + q
            let discriminant = p.w * p.w - 4. * q.w;
            let summand = match discriminant {
                d if d < 0. => quat((-d).sqrt(), 0., 0., 0.),
                d => quat(0., 0., 0., d.sqrt()),
            };

            return ((-p + summand) / 2., (-p - summand) / 2.);
        }

        // q is complex
        let pw2_minus_4qw = p.w * p.w - 4. * q.w;
        let ro = ((pw2_minus_4qw
            + (pw2_minus_4qw * pw2_minus_4qw + 16. * q_im.length_squared()).sqrt())
            / 2.)
            .sqrt();
        let summand = quat(0., 0., 0., ro / 2.) - q_im / ro;
        return (-p / 2. + summand, -p / 2. - summand);
    }

    let b = p_prime.length_squared() + 2. * q_prime.w;
    let e = q_prime.length_squared();
    let d = 2. * p_prime.dot(q_prime);

    let root_fn = |t, n| {
        let p_prime_plus_t = p_prime + quat(0., 0., 0., t);
        quat(0., 0., 0., -p.w / 2.)
            - p_prime_plus_t.conjugate() / p_prime_plus_t.length_squared()
                * (q_prime - quat(0., 0., 0., n))
    };

    let b2_minus_4e = b * b - 4. * e;
    match (b, e, d) {
        (b, _e, d) if d.abs() < EPS && b2_minus_4e >= 0. => {
            let sqrt = b2_minus_4e.sqrt();
            (root_fn(0., (b + sqrt) / 2.), root_fn(0., (b - sqrt) / 2.))
        }
        (b, e, d) if d.abs() < EPS && b2_minus_4e < 0. => {
            let e_sqrt = e.sqrt();
            let t = (2. * e_sqrt - b).sqrt();
            (root_fn(t, e_sqrt), root_fn(-t, e_sqrt))
        }
        (b, _e, d) => {
            let t = roots::find_roots_cubic(1., 2. * b, b2_minus_4e, -d * d)
                .as_ref()
                .iter()
                .copied()
                .find(|&x| x > 0.)
                .expect("polynomial should have a unique positive root")
                .sqrt();

            (
                root_fn(t, (t * t * t + b * t + d) / (2. * t)),
                root_fn(-t, (t * t * t + b * t - d) / (2. * t)),
            )
        }
    }
}

// m-by-N matrix
fn gauss_jordan_eliminate<const N: usize>(rows: &mut [[f32; N]]) {
    const EPS: f32 = 1e-4;
    let m = rows.len();

    let mut h = 0;
    for k in 0..N {
        let i_max = rows[h..]
            .iter()
            .map(|row| OrderedFloat(row[k].abs()))
            .position_max()
            .unwrap()
            + h;
        if rows[i_max][k].abs() < EPS {
            continue;
        }

        rows.swap(h, i_max);
        let (pivot_row, rows_after_pivot) = rows[h..].split_first_mut().unwrap();
        for row in rows_after_pivot {
            let f = row[k] / pivot_row[k];
            row[k] = 0.;
            for j in k + 1..N {
                row[j] -= pivot_row[j] * f;
            }
        }

        h += 1;
        if h >= m {
            break;
        }
    }

    // convert to reduced form via back-substitution
    for i in (0..m).rev() {
        match rows[i].into_iter().find_position(|x| x.abs() > EPS) {
            None => continue,
            Some((k, lead)) => {
                let (pivot_row, rows_before_pivot) = rows[..=i].split_last_mut().unwrap();
                for row in rows_before_pivot {
                    let f = row[k] / pivot_row[k];
                    row[k] = 0.;
                    for j in k + 1..N {
                        row[j] -= pivot_row[j] * f;
                    }
                }

                // normalize the pivot row
                pivot_row[k] = 1.;
                for entry in pivot_row[k + 1..].iter_mut() {
                    *entry /= lead;
                }
            }
        }
    }
}

fn transpose<const M: usize, const N: usize>(mat: [[f32; N]; M]) -> [[f32; M]; N] {
    let mut new = [[0f32; M]; N];
    for (i, row) in new.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = mat[j][i];
        }
    }

    new
}

// Simple polynomial implementation with overflowing multiplication.
// note that N is actually the degree minus one
#[derive(Debug, Clone, Copy)]
struct Polynomial<const N: usize>([f32; N]);

impl<const N: usize> Default for Polynomial<N> {
    fn default() -> Self {
        Self([0f32; N])
    }
}

impl<const N: usize> Polynomial<N> {
    fn identity() -> Self {
        let mut poly = Self::default();
        poly.0[0] = 1.;

        poly
    }

    fn resize<const M: usize>(self) -> Polynomial<M> {
        let mut zeroed = [0f32; M];

        for (i, entry) in zeroed.iter_mut().take(N).enumerate() {
            *entry = self.0[i]
        }

        Polynomial(zeroed)
    }

    fn eval(self, val: f32) -> f32 {
        let mut product = 1.;
        self.0
            .into_iter()
            .map(|coef| {
                let result = coef * product;
                product *= val;
                result
            })
            .sum()
    }
}

impl<const N: usize> Add for Polynomial<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(core::array::from_fn(|k| self.0[k] + rhs.0[k]))
    }
}

impl<const N: usize> AddAssign for Polynomial<N> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const N: usize> Mul for Polynomial<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(core::array::from_fn(|k| {
            let mut sum = 0.;
            for l in 0..k {
                sum += self.0[l] * rhs.0[k - l];
            }
            sum
        }))
    }
}

impl<const N: usize> Mul<f32> for Polynomial<N> {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self(self.0.map(|x| x * rhs))
    }
}

impl<const N: usize> Div<f32> for Polynomial<N> {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self(self.0.map(|x| x / rhs))
    }
}

impl<const N: usize> MulAssign for Polynomial<N> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const N: usize> Shl<usize> for Polynomial<N> {
    type Output = Self;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn shl(mut self, rhs: usize) -> Self::Output {
        for i in rhs..self.0.len() {
            self.0[i] = self.0[i - rhs];
        }

        for i in 0..rhs {
            self.0[i] = 0.;
        }

        self
    }
}

#[cfg(test)]
mod tests {
    use glam::quat;

    use crate::gauss_jordan_eliminate;

    #[test]
    fn quat_quadratic() {
        let a0 = quat(1., 0., 0., 1.);
        let a1 = quat(0., 1., 0., 1.);
        let a2 = quat(0., 0., 1., 2.);

        let (r0, r1) = crate::quaternion_quadratic_roots(a0, a1, a2);
        for root in [r0, r1] {
            let eval = a0 * root * root + a1 * root + a2;
            assert!(eval.length() < 1e-4, "{eval} (eval of {root}) should be 0");
        }
    }

    #[test]
    fn eliminate_seven_by_six() {
        #[allow(clippy::approx_constant)]
        let mut system: [[f32; 7]; 6] = [
            [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, -0.01],
            [1.00, 0.63, 0.39, 0.25, 0.16, 0.10, 0.61],
            [1.00, 1.26, 1.58, 1.98, 2.49, 3.13, 0.91],
            [1.00, 1.88, 3.55, 6.70, 12.62, 23.80, 0.99],
            [1.00, 2.51, 6.32, 15.88, 39.90, 100.28, 0.60],
            [1.00, 3.14, 9.87, 31.01, 97.41, 306.02, 0.02],
        ];

        let solution = [-0.01, 1.60278, -1.61320, 1.24549, -0.49098, 0.06576];

        gauss_jordan_eliminate(&mut system);

        for row in system {
            println!("{row:.2?}");
        }
        let result: [f32; 6] = core::array::from_fn(|i| system[i][6]);
        // assert_eq!();
        for (expected, result) in solution.iter().zip(result.iter()) {
            assert!((expected - result).abs() < 1e-4);
        }
    }
}
