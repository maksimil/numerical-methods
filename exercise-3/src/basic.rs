pub type Scalar = f64;

#[derive(Debug, Clone)]
pub struct IntervalType<S> {
    pub start: S,
    pub end: S,
}

pub type Interval = IntervalType<Scalar>;
pub type IntInterval = IntervalType<i32>;

impl<S> IntervalType<S> {
    pub fn new(start: S, end: S) -> Self {
        Self { start, end }
    }
}

impl IntInterval {
    pub fn diameter(&self) -> u32 {
        (self.end - self.start).unsigned_abs()
    }
}

impl Interval {
    pub fn contains(&self, x: Scalar) -> bool {
        x >= self.start && x <= self.end
    }

    pub fn middle(&self) -> Scalar {
        (self.start + self.end) / Scalar::from(2)
    }
}
