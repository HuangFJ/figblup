use na::*;

/// <pre>
/// A B
/// C D
/// </pre>
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
/// A
/// C
/// </pre>
pub fn compose_v_mat(mut block_a: DMatrix<f64>, mut block_c: DMatrix<f64>) -> DMatrix<f64> {
    assert!(block_a.ncols() == block_c.ncols());

    if block_a.nrows() >= block_c.nrows() {
        let ith = block_a.nrows();
        block_a = block_a.insert_rows(ith, block_c.nrows(), 0.0);

        for i in 0..block_c.nrows() {
            block_a.set_row(ith + i, &block_c.row(i));
        }
        block_a
    } else {
        block_c = block_c.insert_rows(0, block_a.nrows(), 0.0);

        for i in 0..block_a.nrows() {
            block_c.set_row(i, &block_a.row(i));
        }
        block_c
    }
}

/// <pre>
/// A   B
/// </pre>
pub fn compose_h_mat(mut block_a: DMatrix<f64>, block_b: DMatrix<f64>) -> DMatrix<f64> {
    assert!(block_a.nrows() == block_b.nrows());

    block_a.extend(block_b.column_iter());
    block_a
}

/// <pre>
/// a
/// c
/// </pre>
pub fn compose_v_vec(vec_a: DVector<f64>, vec_c: DVector<f64>) -> DVector<f64> {
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

/// <pre>
/// a   b
/// </pre
pub fn compose_h_vec(mut vec_a: RowDVector<f64>, vec_b: RowDVector<f64>) -> RowDVector<f64> {
    vec_a.extend(vec_b.column_iter());
    vec_a
}

#[cfg(all(test))]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let mat_3_3 = dmatrix![
            1.33,2.33,3.33;
            4.33,5.33,6.33;
            7.33,8.33,9.33
        ];
        let mat_3_2 = dmatrix![
            1.32,2.32;
            3.32,4.32;
            5.32,6.32
        ];
        let mat_2_3 = dmatrix![
            1.23,2.23,3.23;
            4.23,5.23,6.23
        ];
        let mat_2_2 = dmatrix![
            1.22,2.22;
            3.22,4.22
        ];

        assert_eq!(
            compose_sqr_mat(
                mat_3_3.clone(),
                mat_3_2.clone(),
                mat_2_3.clone(),
                mat_2_2.clone()
            ),
            dmatrix![
                1.33,2.33,3.33,1.32,2.32;
                4.33,5.33,6.33,3.32,4.32;
                7.33,8.33,9.33,5.32,6.32;
                1.23,2.23,3.23,1.22,2.22;
                4.23,5.23,6.23,3.22,4.22
            ]
        );

        assert_eq!(
            compose_sqr_mat(
                mat_2_2.clone(),
                mat_2_3.clone(),
                mat_3_2.clone(),
                mat_3_3.clone()
            ),
            dmatrix![
                1.22,2.22,1.23,2.23,3.23;
                3.22,4.22,4.23,5.23,6.23;
                1.32,2.32,1.33,2.33,3.33;
                3.32,4.32,4.33,5.33,6.33;
                5.32,6.32,7.33,8.33,9.33
            ]
        );

        assert_eq!(
            compose_v_mat(mat_3_3.clone(), mat_2_3.clone()),
            dmatrix![
                1.33,2.33,3.33;
                4.33,5.33,6.33;
                7.33,8.33,9.33;
                1.23,2.23,3.23;
                4.23,5.23,6.23
            ]
        );

        assert_eq!(
            compose_v_mat(mat_2_2.clone(), mat_3_2.clone()),
            dmatrix![
                1.22,2.22;
                3.22,4.22;
                1.32,2.32;
                3.32,4.32;
                5.32,6.32
            ]
        );

        assert_eq!(
            compose_h_mat(mat_2_3.clone(), mat_2_2.clone()),
            dmatrix![
                1.23,2.23,3.23,1.22,2.22;
                4.23,5.23,6.23,3.22,4.22
            ]
        );

        assert_eq!(
            compose_h_vec(mat_3_3.row(0).clone_owned(), mat_3_3.row(1).clone_owned()),
            dmatrix![1.33, 2.33, 3.33, 4.33, 5.33, 6.33]
        );

        assert_eq!(
            compose_v_vec(
                mat_2_2.column(0).into_owned(),
                mat_3_3.column(0).into_owned()
            ),
            dmatrix![
                1.22;
                3.22;
                1.33;
                4.33;
                7.33
            ]
        );
    }
}
