pub type Index = usize;
pub type Scalar = f64;

pub const SCALAR_ZERO: Scalar = 0.;
pub const SCALAR_ONE: Scalar = 1.;

pub const MAX_ITERATIONS: Index = 1_000_000;

// LU decomposition
pub const MIN_PIVOT: Scalar = 1e-7;

// Random matrix generation
pub const VALUES_SEPARATION: Scalar = 1e-1;
pub const VALUES_MAGNITUDE_MUL: Scalar = 100.; // should be > 1
pub const RANDOM_TRANSFROMATION_VALUES_MAGNITUDE: Scalar = 10.;

// Power method
pub const PM_MIN_COORD: Scalar = 1e-8;
pub const PM_CONVERGENCE: Scalar = 1e-6;

// QR-algorithm
pub const SUBDIAGONAL_SMALL: Scalar = 1e-8;
