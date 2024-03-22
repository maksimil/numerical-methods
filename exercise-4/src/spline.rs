use std::cmp::Ordering;

use crate::{lu::LUDecomposition, polynomial::poly_at, scalar::*};

#[derive(Debug, Clone)]
pub struct Spline {
    coefs: Vec<Scalar>,

    points: Vec<Scalar>,
    order: Index,
}

impl Spline {
    pub fn at(&self, x: Scalar) -> Scalar {
        let idx = {
            let idx = self
                .points
                .binary_search_by(|y| {
                    if y < &x {
                        Ordering::Less
                    } else {
                        Ordering::Greater
                    }
                })
                .err()
                .unwrap();
            if idx == 0 {
                1
            } else if idx == self.points.len() {
                self.points.len() - 1
            } else {
                idx
            }
        };

        poly_at(
            &self.coefs[(idx - 1) * (self.order + 1)..idx * (self.order + 1)],
            x,
        )
    }

    // points array should be sorted
    pub fn new(
        order: Index,
        points: &Vec<Scalar>,
        f: impl Fn(Scalar) -> Scalar,
        initial: &Vec<Scalar>, // derivatives 1..=order-1, len=order-1
    ) -> Spline {
        // computing the needed factorials
        let factorial = {
            let factorial_to = order;
            let mut factorial = vec![SCALAR_ZERO; factorial_to + 1];
            factorial[0] = SCALAR_ONE;
            for n in 1..=factorial_to {
                factorial[n] = (n as Scalar) * factorial[n - 1];
            }
            factorial
        };

        let k = points.len() - 1;

        let rowlen = (order + 1) * k;
        let mut coefs = vec![SCALAR_ZERO; rowlen]; // a_{jn} = coefs[j*order+n]
        let mut matrix = vec![SCALAR_ZERO; rowlen * rowlen];

        // filling the linear system
        {
            // values at points (rows 0..2k)
            for j in 0..k {
                // l_j(x_j)=y_j
                let mut xjn = SCALAR_ONE;
                for n in 0..=order {
                    matrix[(2 * j) * rowlen + (order + 1) * j + n] = xjn;
                    xjn *= points[j];
                }
                coefs[2 * j] = f(points[j]);

                // lj(x_{j+1})=y_{j+1}
                let mut xjn = SCALAR_ONE;
                for n in 0..=order {
                    matrix[(2 * j + 1) * rowlen + (order + 1) * j + n] = xjn;
                    xjn *= points[j + 1];
                }
                coefs[2 * j + 1] = f(points[j + 1]);
            }

            // initial conditions (rows 2k..2k+order-1)
            for i in 1..=(order - 1) {
                // 0..i are zero
                // i is i!
                matrix[(2 * k + i - 1) * rowlen + i] = factorial[i];
                // i+1..=order is l!/(l-i)!*x^{l-i}
                let mut xjn = SCALAR_ONE;
                for n in (i + 1)..=order {
                    xjn *= points[0];
                    matrix[(2 * k + i - 1) * rowlen + n] = xjn * factorial[n] / factorial[n - i];
                }
                coefs[2 * k + i - 1] = initial[i - 1];
            }

            // continuous derivative (rows 2k+order-1..2k+k*(order+1))
            for i in 1..=(order - 1) {
                for j in 0..k - 1 {
                    // 0..i are zero
                    // i is i!
                    matrix[((2 * k + order - 1) + j * (order - 1) + i - 1) * rowlen
                        + j * (order + 1)
                        + i] = factorial[i];
                    matrix[((2 * k + order - 1) + j * (order - 1) + i - 1) * rowlen
                        + (j + 1) * (order + 1)
                        + i] = -factorial[i];
                    // i+1=order is l!/(l-i)!*x^{l-i}
                    let mut xjn = SCALAR_ONE;
                    for n in (i + 1)..=order {
                        xjn *= points[j + 1];
                        matrix[((2 * k + order - 1) + j * (order - 1) + i - 1) * rowlen
                            + j * (order + 1)
                            + n] = xjn * factorial[n] / factorial[n - i];
                        matrix[((2 * k + order - 1) + j * (order - 1) + i - 1) * rowlen
                            + (j + 1) * (order + 1)
                            + n] = -xjn * factorial[n] / factorial[n - i];
                    }
                }
            }
        }

        println!("dim={rowlen}");

        // computing coefs
        let lud = LUDecomposition::compute(&matrix, rowlen);
        lud.solve(&mut coefs);

        Spline {
            coefs,
            points: points.clone(),
            order,
        }
    }
}
