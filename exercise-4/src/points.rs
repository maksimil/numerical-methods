use crate::scalar::*;

pub fn points_eq(a: Scalar, b: Scalar, n: Index) -> Vec<Scalar> {
    let mut c = vec![SCALAR_ZERO; n + 1];

    for k in 0..=n {
        c[k] = a + (b - a) * (k as Scalar) / (n as Scalar)
    }

    c
}

pub fn points_cheb(a: Scalar, b: Scalar, n: Index) -> Vec<Scalar> {
    let mut c = vec![SCALAR_ZERO; n + 1];

    for k in 0..=n {
        c[k] = 0.5
            * ((b - a)
                * (SCALAR_PI * (2. * ((n - k) as Scalar) + 1.) / (2. * ((n as Scalar) + 1.)))
                    .cos()
                + (b + a));
    }

    c
}
