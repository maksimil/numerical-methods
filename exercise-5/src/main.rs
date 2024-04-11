use rand_matrix::generate_matrix;

mod lu;
mod rand_matrix;
mod scalar;

fn main() {
    let dimension = 10;
    let generated = generate_matrix(dimension);

    println!("{:?}", generated);
}
