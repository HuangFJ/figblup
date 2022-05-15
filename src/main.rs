extern crate nalgebra as na;

mod nrm;
// use std::intrinsics::sqrtf64;

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

fn matA(inds: Vec<(i64, i64, i64)>) -> DMatrix<f64> {
    /// A = TDT'
    /// D是孟德尔抽样的方差协方差矩阵，个体的育种值 u_i，亲本的育种值是u_s,u_d，则孟德尔抽样 m_i=u_i-0.5*(u_s+u_d)
    /// d_ii = var(m_i)/sigma_u^2 = 0.5-0.25*(Fs+Fd)
    /// 双亲未知 d_ii = 1
    ///
    /// T的t_ij是个体i与个体j的亲缘系数，并假设不存在近交，从一个世代到下一个世代追溯基因的流动，只需解释亲本到后代的关系
    /// t_ii = 1; t_ij=0.5*(t_sj+t_dj); t_ij=0.5*(t_sj); t_ij=0
    ///
    /// A对角元素 a_ii = 1+Fi  (Fi是个体i得近交系数)
    /// 如果u_i 是个体育种值，则var(u_i) = a_ii * sigma_u^2=(1+Fi)*sigma_u^2
    ///
    ///
    let n = inds.len();
    let mut A = DMatrix::<f64>::zeros(n, n);
    for ind in inds {
        let ith = (ind.0 - 1) as usize;
        if ind.1 > 0 && ind.2 > 0 {
            //双亲都已知
            for jth in 0..ith {
                A[(ith, jth)] =
                    0.5 * (A[(jth, (ind.1 - 1) as usize)] + A[(jth, (ind.2 - 1) as usize)]);
                A[(jth, ith)] = A[(ith, jth)];
            }
            A[(ith, ith)] = 1.0 + 0.5 * A[((ind.1 - 1) as usize, (ind.2 - 1) as usize)];
        } else if ind.1 == 0 && ind.2 == 0 {
            //双亲均未知，且双亲不相关
            for jth in 0..ith {
                A[(ith, jth)] = 0.0;
                A[(jth, ith)] = A[(ith, jth)];
            }
            A[(ith, ith)] = 1.0;
        } else {
            //仅知一个亲本，并假设与配偶不相关
            let p = if ind.1 > 0 { ind.1 } else { ind.2 };
            for jth in 0..ith {
                A[(ith, jth)] = 0.5 * A[(jth, (p - 1) as usize)];
                A[(jth, ith)] = A[(ith, jth)];
            }
            A[(ith, ith)] = 1.0;
        }
    }
    A
}

fn matDI(inds: Vec<(i64, i64, i64)>) -> DMatrix<f64>{
    /// L=T*sqrt(D)
    /// A=LL'
    /// l_ii=sqrt(1.0-0.25*(sum_m_1-s(l_sm^2)+sum_m_1-d(l_dm^2)))
    /// DI的对角线元素a_i，a_i=1/l_ii^2
    let n = inds.len();
    let mut DI = DMatrix::<f64>::zeros(n, n);
    let mut L = DMatrix::<f64>::zeros(n, n);
    L[(0, 0)] = 1.0;
    for ind in inds {
        let ith = (ind.0 - 1) as usize;
        for jth in 0..ith {
            let s = if jth as i64 <= (ind.1 - 1) {
                L[((ind.1 - 1) as usize, jth)]
            } else {
                0.0
            };
            let d = if jth as i64 <= (ind.2 - 1) {
                L[((ind.2 - 1) as usize, jth)]
            } else {
                0.0
            };
            L[(ith, jth)] = 0.5 * (s + d);
        }

        let sum_s = (0..(ind.1))
            .map(|s_idx| {
                let val = L[((ind.1 - 1) as usize, s_idx as usize)];
                val * val
            })
            .reduce(|accum, item| accum + item)
            .unwrap_or(0.0);
        let sum_d = (0..(ind.2))
            .map(|d_idx| {
                let val = L[((ind.2 - 1) as usize, d_idx as usize)];
                val * val
            })
            .reduce(|accum, item| accum + item)
            .unwrap_or(0.0);

        L[(ith, ith)] = (1.0 - 0.25 * (sum_s + sum_d)).sqrt();
        
        DI[(ith, ith)] = 1.0 / (L[(ith, ith)]).powi(2);
    }

    println!("{}", L);
    DI
}

fn matAI(inds: Vec<(i64, i64, i64)>, matDI: DMatrix<f64>) -> DMatrix<f64> {
    ///
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
            //双亲都已知，没有近交的DI的对角线元素
            let diag = matDI[diag_idx];
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

    let pedi = vec![
        (1, 0, 0),
        (2, 0, 0),
        (3, 1, 2),
        (4, 1, 0),
        (5, 4, 3),
        (6, 5, 2),
    ];
    let DI = matDI(pedi);

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
    let mat = matAI(pedi, DI);
    println!("{}", &mat);

    // let pedi = vec![
    //     (1, 0, 0),
    //     (2, 0, 0),
    //     (3, 1, 2),
    //     (4, 1, 0),
    //     (5, 4, 3),
    //     (6, 5, 2),
    // ];
    // let mat = matA(pedi);
    // println!("{}", &mat);
}
