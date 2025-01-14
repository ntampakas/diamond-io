/// Print a vector of integers
pub fn print_vector(label: &str, vec: &Vec<u64>) {
    print!("{} [", label);
    for (k, &val) in vec.iter().enumerate() {
        if k > 0 {
            print!(", ");
        }
        print!("{}", val);
    }
    println!("]");
}

/// Print a matrix of ring elements
pub fn print_matrix_ring(label: &str, matrix: &Vec<Vec<Vec<u64>>>) {
    println!("\n{}", label,);

    for (i, row) in matrix.iter().enumerate() {
        for (j, col) in row.iter().enumerate() {
            print!("r{}c{}: ", i, j);
            print_vector("", col);
        }
    }
}

/// Print a vector of ring elements
pub fn print_vector_ring(label: &str, vec: &Vec<Vec<u64>>) {
    println!("\n{}", label);
    for (i, inner_vec) in vec.iter().enumerate() {
        print!("{}[{}]: [", label, i);
        for (j, &val) in inner_vec.iter().enumerate() {
            if j > 0 {
                print!(", ");
            }
            print!("{}", val);
        }
        println!("]");
    }
}
