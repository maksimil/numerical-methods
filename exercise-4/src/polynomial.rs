use std::num::Wrapping;

pub type Scalar = f64;
pub const SCALAR_ZERO: Scalar = 0.;
pub const SCALAR_ONE: Scalar = 1.;

pub type Index = usize;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coefs: Vec<Scalar>,
}

impl Polynomial {
    pub fn new(n: Index) -> Polynomial {
        let coefs = vec![SCALAR_ZERO; n + 1];
        Polynomial { coefs }
    }

    pub fn at(&self, x: Scalar) -> Scalar {
        let mut y = SCALAR_ZERO;
        let mut xn = SCALAR_ONE;

        for n in 0..self.coefs.len() {
            y += xn * self.coefs[n];
            xn *= x;
        }

        y
    }
}

pub fn add_root(root: Scalar, mut p: Polynomial) -> Polynomial {
    let mut store = SCALAR_ZERO;

    for n in 0..p.coefs.len() {
        let c_prev = store;
        store = p.coefs[n];
        p.coefs[n] = c_prev - root * p.coefs[n];
    }

    p
}

fn create_lk(points: &Vec<Scalar>, index: Index) -> Polynomial {
    let mut p = Polynomial::new(points.len() - 1);
    p.coefs[0] = SCALAR_ONE;

    for k in 0..points.len() {
        if k == index {
            continue;
        }

        p = add_root(points[k], p);
    }

    let d = p.at(points[index]);

    for n in 0..p.coefs.len() {
        p.coefs[n] /= d;
    }

    p
}

pub fn create_lagrange(points: &Vec<Scalar>, f: impl Fn(Scalar) -> Scalar) -> Polynomial {
    let mut p = Polynomial::new(points.len() - 1);

    for k in 0..points.len() {
        let value = f(points[k]);
        let lk = create_lk(points, k);
        for n in 0..p.coefs.len() {
            p.coefs[n] += lk.coefs[n] * value
        }
    }

    p
}

fn diff_table_index(n: Index, k: Index) -> Index {
    k * n - k * (Wrapping(k) - Wrapping(1)).0 / 2
}

pub fn create_newton(points: &Vec<Scalar>, f: impl Fn(Scalar) -> Scalar) -> Polynomial {
    // fill the differences table
    let n = points.len();
    let mut diff_table = vec![SCALAR_ZERO; diff_table_index(n, n - 1) + 1];

    for k in 0..n {
        diff_table[k] = f(points[k]);
    }

    for k in 1..n {
        let ck_prev = diff_table_index(n, k - 1);
        let ck = diff_table_index(n, k);

        for i in 0..n - k {
            diff_table[ck + i] = (diff_table[ck_prev + i + 1] - diff_table[ck_prev + i])
                / (points[i + k] - points[i]);
        }
    }

    // construct the polynomial
    let mut p = Polynomial::new(n - 1);
    let mut omega = Polynomial::new(n - 1);
    omega.coefs[0] = SCALAR_ONE;

    for k in 0..n {
        let value = diff_table[diff_table_index(n, k)];
        for i in 0..n {
            p.coefs[i] += value * omega.coefs[i];
        }

        omega = add_root(points[k], omega);
    }

    p
}
