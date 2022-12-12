# pythagorean-hodographs

[![Build Status][actions-badge]][actions-url]
![License][license-badge]
[![Crates.io][crates-badge]][crates-url]
[![API][docs-badge]][docs-url]

[actions-badge]: https://github.com/suremarc/pythagorean-hodographs/workflows/build/badge.svg?event=push
[actions-url]: https://github.com/suremarc/pythagorean-hodographs/actions?query=workflow%3Abuild+branch%3Amaster
[docs-badge]: https://docs.rs/pythagorean-hodographs/badge.svg
[docs-url]: https://docs.rs/pythagorean-hodographs
[license-badge]: https://img.shields.io/badge/license-MIT_OR_Apache--2.0-blue.svg
[crates-badge]: https://img.shields.io/crates/v/pythagorean-hodographs.svg
[crates-url]: https://crates.io/crates/pythagorean-hodographs

Implementations of Pythagorean hodograph (PH) splines in Rust, based on the [work of R. Farouki](https://doi.org/10.1007/978-3-540-73398-0) and [this paper by Zbyněk Šír & Bert Jüttler](https://doi.org/10.1007/11537908_22).

## What are Pythagorean hodographs?

Pythaogrean hodographs are interpolating polynomial curves that satisfy an analogue of the Pythagorean identity. To illustrate, let $\mathbf{p}(t)$ be a polynomial parametric curve in $\mathbb{R}^2$. If there exist polynomials $u(t), v(t), w(t)$ such that the following holds:

$$\begin{equation}
    \frac{d\mathbf{p}(t)}{dt} = \begin{bmatrix}
        u(t)^2 - v(t)^2 \\
        2u(t)v(t) \\
    \end{bmatrix}\cdot w(t)
\end{equation}$$

then $\mathbf{p}(t)$ is a Pythaogrean hodograph, because $(u^2-v^2)^2+(2uv)^2=(u^2+v^2)^2,$ a perfect square of a polynomial. It then follows that $|\mathbf{p}'(t)|=(u(t)^2+v(t)^2)\cdot w(t).$ $\mathbf{p}(t)$ is then a __planar PH curve.__

There are many useful properties of Pythagorean hodographs as a result of this:

1. The tangent vector, $\frac{\mathbf{p}'(t)}{|\mathbf{p}'(t)|},$ has rational functions as components. This makes computing the tangent fast, although square roots nowadays have dedicated CPU instructions.

2. The arc-length function $s(t)=\int_0^t|\mathbf{p}'(\tau)|\ d\tau$ is a polynomial function.

3. As a corollary of #2, computing the arc-length parameterization is trivial using Newton's method.

4. In three dimensions, Pythagorean hodographs have a naturally associated moving frame called the _Euler-Rodrigues frame._

## Spatial PH curves

Well-studied readers may have noticed a connection to complex numbers. In particular, if $x(t)=u(t)+iv(t),$ then $p(t)=x(t)^2$ (where $\mathbf{p}(t)$ and $p(t)$ are equivalent via the map from $\mathbb{R}^2\rightarrow\mathbb{C}$).

Complex numbers are isomorphic to the plane, but we need something else if we want to model three dimensions -- __quaternions__. One has the following presentation: $\mathbf{r}(t)$ is a spatial PH curve if there exists a quaternion curve $\mathcal{A}(t)$ such that

$$\begin{equation}
\begin{gathered}
    \begin{aligned}
    \mathcal{A}(t)&=u(t)+v(t)\cdot\mathbf{i}+p(t)\cdot\mathbf{j}+q(t)\cdot\mathbf{k} \\
    r'(t)&=\mathcal{A}(t)\ \mathbf{i}\ \bar{\mathcal{A}}(t)\\
    &=[u(t)^2+v(t)^2+p(t)^2+q(t)^2]\cdot\mathbf{i}+2[u(t)q(t)+v(t)p(t)]\cdot\mathbf{j}+2[v(t)q(t)-u(t)p(t)]\cdot\mathbf{k}
    \end{aligned}
\end{gathered}
\end{equation}$$

It then holds that $|\mathbf{r}'(t)|=|\mathcal{A}(t)|^2.$

### The Euler-Rodrigues frame

With spatial PH curves, there is a naturally associated moving frame called the Euler-Rodrigues frame. One has the following formula for the tangent, normal, and binormal:

$$\begin{equation}
    \begin{bmatrix}
    \mathbf{T}(t) \\
    \mathbf{N}(t) \\
    \mathbf{B}(t)
    \end{bmatrix} = \begin{bmatrix}
    \mathcal{A}(t)\ \mathbf{i}\ \bar{\mathcal{A}}(t)\\
    \mathcal{A}(t)\ \mathbf{j}\ \bar{\mathcal{A}}(t) \\
    \mathcal{A}(t)\ \mathbf{k}\ \bar{\mathcal{A}}(t)
    \end{bmatrix}
\end{equation}$$

or, more succinctly, the rotation matrix $[\mathbf{T\ N\ B}]$ is equivalent to the 3-dimensional rotation associated with the quaternion $\mathcal{A}(t).$
