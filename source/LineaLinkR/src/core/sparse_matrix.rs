#[derive(Debug, Clone)]
pub struct SparseMatrix {
    pub rows: Vec<usize>,
    pub col_ptrs: Vec<usize>,
    pub values: Vec<f32>,
    pub nrows: usize,
    pub ncols: usize,
}

impl SparseMatrix {
    pub fn new_csc(
        rows: Vec<usize>,
        col_ptrs: Vec<usize>,
        values: Vec<f32>,
        nrows: usize,
        ncols: usize,
    ) -> Self {
        Self {
            rows,
            col_ptrs,
            values,
            nrows,
            ncols,
        }
    }
}
