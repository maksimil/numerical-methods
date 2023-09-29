fn taylor_sin(acc: f64, x: f64) -> f64 {
    let mut uk = x;
    let mut r = uk;
    let mut k = 1.;

    while uk.abs() > acc {
        uk *= -(x * x) / ((2. * k + 1.) * (2. * k));
        k += 1.;
        r += uk;
    }

    r
}

fn geron_sqrt(acc: f64, x: f64) -> f64 {
    let mut pthis = x + 0.5;
    let mut pnext = (pthis + (x / pthis)) / 2.;

    while (pnext - pthis).abs() > acc {
        pthis = pnext;
        pnext = (pnext + (x / pnext)) / 2.;
    }

    pnext
}

struct Computation {
    u: f64,
    v: f64,
    z: f64,
}

fn approx_z(acc: f64, x: f64) -> Computation {
    let sin_acc = acc / 2.5;
    let sqrt_acc = acc / 2.;

    let u = taylor_sin(sin_acc, 4.5 * x + 0.6);
    let v = geron_sqrt(sqrt_acc, 1.0 + x - 12. * x * x);
    let z = u / v;

    Computation { u, v, z }
}

fn compare_z(x: f64) -> Computation {
    let u = (4.5 * x + 0.6).sin();
    let v = (1.0 + x - 12. * x * x).sqrt();
    let z = u / v;

    Computation { u, v, z }
}

const EPSILON: f64 = 0.00_000_1;
const U_ACC: f64 = EPSILON / 2.5;
const V_ACC: f64 = EPSILON / 2.;

fn main() {
    let x_start = 0.1;
    let x_end = 0.2;
    let x_step = 0.01;
    let n_max = (((x_end - x_start) / x_step) as f64).ceil() as usize;

    println!(
        "x,u,\\Delta_u,\\bar{{u}},\\bar{{\\Delta_u}},\
             v,\\Delta_v,\\bar{{v}},\\bar{{\\Delta_v}},\
             z,\\Delta_z,\\bar{{z}},\\bar{{\\Delta_z}}"
    );

    for n in 0..=n_max {
        let x = x_start + (n as f64) * x_step;
        let approx_computation = approx_z(EPSILON, x);
        let compare_computation = compare_z(x);

        let u = approx_computation.u;
        let v = approx_computation.v;
        let z = approx_computation.z;

        let bu = compare_computation.u;
        let bv = compare_computation.v;
        let bz = compare_computation.z;

        let du = U_ACC;
        let dv = V_ACC;
        let dz = EPSILON;

        let bdu = bu - u;
        let bdv = bv - v;
        let bdz = bz - z;

        println!(
            "{x},{u},{du},{bu},{bdu},\
                 {v},{dv},{bv},{bdv},\
                 {z},{dz},{bz},{bdz}"
        );
    }
}
