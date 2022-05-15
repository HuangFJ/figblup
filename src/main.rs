extern crate nalgebra as na;

mod nrm;
use na::*;

fn consist4matrixs(
    mut matA: DMatrix<f64>,
    mut matB: DMatrix<f64>,
    mut matC: DMatrix<f64>,
    mut matD: DMatrix<f64>,
) -> DMatrix<f64> {
    /// A B
    /// C D
    // let R = consist4matrixs(
    //     dmatrix![
    //         1.1,2.0,3.0;
    //         4.0,5.0,6.0;
    //         7.0,8.0,9.1
    //     ],
    //     dmatrix![
    //         1.2,2.0;
    //         4.0,5.0;
    //         7.0,8.2
    //     ],
    //     dmatrix![
    //         1.3,2.0,3.0;
    //         4.0,5.0,6.3
    //     ],
    //     dmatrix![
    //         1.4,2.0;
    //         4.0,5.4
    //     ],
    // );
    // println!("{}", R);

    // let R = consist4matrixs(
    //     dmatrix![
    //         1.2,2.0;
    //         4.0,5.2
    //     ],
    //     dmatrix![
    //         1.1,2.0,3.0;
    //         4.0,5.0,6.1
    //     ],
    //     dmatrix![
    //         1.4,2.0;
    //         4.0,5.0;
    //         7.0,8.4
    //     ],
    //     dmatrix![
    //         1.3,2.0,3.0;
    //         4.0,5.0,6.0;
    //         7.0,8.0,9.3
    //     ],

    // );
    assert!(matA.nrows() == matB.nrows());
    assert!(matC.nrows() == matD.nrows());
    assert!(matA.ncols() == matC.ncols());
    assert!(matB.ncols() == matD.ncols());

    let (mut matL, mut matR) = if matA.nrows() >= matC.nrows() {
        let ith = matA.nrows();
        matA = matA.insert_rows(ith, matC.nrows(), 0.0);
        matB = matB.insert_rows(ith, matC.nrows(), 0.0);

        for i in 0..matC.nrows() {
            matA.set_row(ith + i, &matC.row(i));
            matB.set_row(ith + i, &matD.row(i));
        }
        (matA, matB)
    } else {
        matC = matC.insert_rows(0, matA.nrows(), 0.0);
        matD = matD.insert_rows(0, matA.nrows(), 0.0);

        for i in 0..matA.nrows() {
            matC.set_row(i, &matA.row(i));
            matD.set_row(i, &matB.row(i));
        }
        (matC, matD)
    };

    //matL.extend(matR.column_iter());

    if matL.ncols() >= matR.ncols() {
        let ith = matL.ncols();
        matL = matL.insert_columns(ith, matR.ncols(), 0.0);
        for i in 0..matR.ncols() {
            matL.set_column(ith + i, &matR.column(i));
        }
        matL
    } else {
        matR = matR.insert_columns(0, matL.ncols(), 0.0);
        for i in 0..matL.ncols() {
            matR.set_column(i, &matL.column(i));
        }
        matR
    }
}

fn matAI_without_jj(inds: Vec<(i64, i64, i64)>) -> DMatrix<f64> {
    /// (1, 0, 0)
    /// (2, 0, 0)
    /// (3, 1, 2)
    /// (4, 1, 0)
    /// (5, 4, 3)
    /// (6, 5, 2)
    ///
    ///
    ///
    let n = inds.len();
    let mut AI = DMatrix::<f64>::zeros(n, n);
    for ind in inds {
        let diag_idx = ((ind.0 - 1) as usize, (ind.0 - 1) as usize);
        if ind.1 > 0 && ind.2 > 0 {
            //双亲都已知
            let diag = 2.0;
            AI[diag_idx] += diag;

            //亲子加性
            let c_p_addtioner = -diag / 2.0;
            for c_p_idx in [
                ((ind.0 - 1) as usize, (ind.1 - 1) as usize),
                ((ind.0 - 1) as usize, (ind.2 - 1) as usize),
                ((ind.1 - 1) as usize, (ind.0 - 1) as usize),
                ((ind.2 - 1) as usize, (ind.0 - 1) as usize),
            ] {
                AI[c_p_idx] += c_p_addtioner;
            }

            //双亲加性
            let p_p_addtioner = diag / 4.0;
            for p_p_idx in [
                ((ind.1 - 1) as usize, (ind.1 - 1) as usize),
                ((ind.2 - 1) as usize, (ind.2 - 1) as usize),
                ((ind.1 - 1) as usize, (ind.2 - 1) as usize),
                ((ind.2 - 1) as usize, (ind.1 - 1) as usize),
            ] {
                AI[p_p_idx] += p_p_addtioner;
            }
        } else if ind.1 == 0 && ind.2 == 0 {
            //双亲均未知
            let diag = 1.0;
            AI[diag_idx] += diag;
        } else {
            //仅知一个亲本
            let diag = 4.0 / 3.0;
            AI[diag_idx] += diag;

            let p = if ind.1 > 0 { ind.1 } else { ind.2 };
            //亲子加性
            let c_p_addtioner = -diag / 2.0;
            for c_p_idx in [
                ((ind.0 - 1) as usize, (p - 1) as usize),
                ((p - 1) as usize, (ind.0 - 1) as usize),
            ] {
                AI[c_p_idx] += c_p_addtioner;
            }

            //双亲加性
            let p_p_addtioner = diag / 4.0;
            let p_p_idx = ((p - 1) as usize, (p - 1) as usize);
            AI[p_p_idx] += p_p_addtioner;
        }
    }
    AI
}

fn main() {
    //| XT*X  XT*Z            ||b| = |XT*Y|
    //| ZT*X  ZT*Z+AI*alpha   ||a| = |ZT*Y|
    //假设 sigma_a_power=20, sigma_e_power=40
    //则 alpha=sigma_e_power/sigma_a_power=2
    //Yij=Pi+Aj+Eij
    //Yij=i性别j个体的WWG；Pi=i性别的固定效应；Aj=j个体的随机效应；Eij=随机误差效应
    //
    //ID    Sex BoarID    SowID WWG（kg）
    //4 公  1   Nil 4.5
    //5 母  3   2   2.9
    //6 母  1   2   3.9
    //7 公  4   5   3.5
    //8 公  3   6   5.0

    //固定效应（性别）矩阵X
    //  公    母
    //4 1   0
    //5 0   1
    //6 0   1
    //7 1   0
    //8 1   0
    let X = dmatrix![
        1.0,  0.0;
        0.0,  1.0;
        0.0,  1.0;
        1.0,  0.0;
        1.0,  0.0;
    ];

    //Z是记录个体与所有个体的数据有无情况，不管有没有数据。1-3的个体是亲本，没有记录
    /*
        1   2   3   4   5   6   7   8
    4               1
    5                   1
    6                       1
    7                           1
    8                               1
    */
    let Z = dmatrix![
        0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0;
    ];
    let Y = dvector![4.5, 2.9, 3.9, 3.5, 5.0];

    let XT_X = &X.transpose() * &X;
    let XT_Z = &X.transpose() * &Z;
    let ZT_X = XT_Z.transpose();
    let ZT_Z = &Z.transpose() * &Z;
    let XT_Y = &X.transpose() * &Y;
    let ZT_Y = &Z.transpose() * &Y;

    let R = consist4matrixs(XT_X, XT_Z, ZT_X, ZT_Z);

    println!("{}", R);

    // println!("{}", XT_X);

    // let diag = dvector![40.0, 40.0, 40.0, 50.0, 50.0];
    // let r = DMatrix::<f32>::from_diagonal(&diag);
    // let r_inverse = r.try_inverse().unwrap();

    // let w = dmatrix![
    //     0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    //     0.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    //     0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    //     0.0, 0.0, 0.0, 0.5, 0.5, 0.0;
    //     0.0, 0.0, 0.5, 0.0, 0.0, 0.5;
    // ];

    // let wT = w.transpose();
    // let result = wT * r_inverse * w;
    // println!("{}", result);
    let pedi = vec![
        (1, 0, 0),
        (2, 0, 0),
        (3, 1, 2),
        (4, 1, 0),
        (5, 4, 3),
        (6, 5, 2),
    ];
    let mat = matAI_without_jj(pedi);
    println!("{}", &mat);
}
