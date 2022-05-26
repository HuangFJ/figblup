use na::*;

/// <pre>
/// A B
/// C D
/// </pre>
/// ```
/// let r = compose_sqr_mat(
///     dmatrix![
///         1.1,2.0,3.0;
///         4.0,5.0,6.0;
///         7.0,8.0,9.1
///     ],
///     dmatrix![
///         1.2,2.0;
///         4.0,5.0;
///         7.0,8.2
///     ],
///     dmatrix![
///         1.3,2.0,3.0;
///         4.0,5.0,6.3
///     ],
///     dmatrix![
///         1.4,2.0;
///         4.0,5.4
///     ],
/// );
/// println!("{}", R);
/// ```
/// ```
/// let r = compose_sqr_mat(
///     dmatrix![
///         1.2,2.0;
///         4.0,5.2
///     ],
///     dmatrix![
///         1.1,2.0,3.0;
///         4.0,5.0,6.1
///     ],
///     dmatrix![
///         1.4,2.0;
///         4.0,5.0;
///         7.0,8.4
///     ],
///     dmatrix![
///         1.3,2.0,3.0;
///         4.0,5.0,6.0;
///         7.0,8.0,9.3
///     ],
/// );
/// ```
pub fn compose_sqr_mat(
    mut block_a: DMatrix<f64>,
    mut block_b: DMatrix<f64>,
    mut block_c: DMatrix<f64>,
    mut block_d: DMatrix<f64>,
) -> DMatrix<f64> {
    assert!(block_a.nrows() == block_b.nrows());
    assert!(block_c.nrows() == block_d.nrows());
    assert!(block_a.ncols() == block_c.ncols());
    assert!(block_b.ncols() == block_d.ncols());

    let (block_l, block_r) = if block_a.nrows() >= block_c.nrows() {
        let ith = block_a.nrows();
        block_a = block_a.insert_rows(ith, block_c.nrows(), 0.0);
        block_b = block_b.insert_rows(ith, block_c.nrows(), 0.0);

        for i in 0..block_c.nrows() {
            block_a.set_row(ith + i, &block_c.row(i));
            block_b.set_row(ith + i, &block_d.row(i));
        }
        (block_a, block_b)
    } else {
        block_c = block_c.insert_rows(0, block_a.nrows(), 0.0);
        block_d = block_d.insert_rows(0, block_a.nrows(), 0.0);

        for i in 0..block_a.nrows() {
            block_c.set_row(i, &block_a.row(i));
            block_d.set_row(i, &block_b.row(i));
        }
        (block_c, block_d)
    };

    //block_l.extend(block_r.column_iter());

    if block_l.ncols() >= block_r.ncols() {
        let ith = block_l.ncols();
        let mut mat = block_l.insert_columns(ith, block_r.ncols(), 0.0);
        for i in 0..block_r.ncols() {
            mat.set_column(ith + i, &block_r.column(i));
        }
        mat
    } else {
        let mut mat = block_r.insert_columns(0, block_l.ncols(), 0.0);
        for i in 0..block_l.ncols() {
            mat.set_column(i, &block_l.column(i));
        }
        mat
    }
}

/// <pre>
/// a
/// c
/// </pre>
pub fn compose_v_vec(mut vec_a: DVector<f64>, mut vec_c: DVector<f64>) -> DVector<f64> {
    if vec_a.nrows() >= vec_c.nrows() {
        let ith = vec_a.nrows();
        let mut vec = vec_a.insert_rows(ith, vec_c.nrows(), 0.0);

        for i in 0..vec_c.nrows() {
            vec.set_row(ith + i, &vec_c.row(i));
        }
        vec
    } else {
        let mut vec = vec_c.insert_rows(0, vec_a.nrows(), 0.0);

        for i in 0..vec_a.nrows() {
            vec.set_row(i, &vec_a.row(i));
        }
        vec
    }
}
