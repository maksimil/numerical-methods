use crate::{lu::LUDecomposition, scalar::*};

#[derive(Debug, Clone)]
pub struct Spline {
    coefs: Vec<Scalar>,

    points: Vec<Scalar>,
    order: Index,
}

impl Spline {
    pub fn new(order: Index, points: Vec<Scalar>, f: impl Fn(Scalar) -> Scalar) -> Spline {
        let k = points.len() - 1;

        let rowlen = (order + 1) * k;
        let mut coefs = vec![SCALAR_ZERO; rowlen]; // a_{jn} = coefs[j*order+n]
        let mut matrix = vec![SCALAR_ZERO; rowlen * rowlen];

        // filling the linear system
        {
            // values at points
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

            assert!(order == 1);
            // initial conditions
            // continuous derivative
        }

        println!("m={matrix:?}");
        println!("y={coefs:?}");

        // computing coefs
        let lud = LUDecomposition::compute(&matrix, rowlen);
        lud.solve(&mut coefs);

        Spline {
            coefs,
            points,
            order,
        }
    }
}
