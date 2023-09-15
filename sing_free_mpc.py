from casadi import *
from acados_template import *
from constants import *
import scipy.linalg
from icecream import ic

__all__ = ["gen_solver"]


def gen_solver(
    model_name: str, s_ref: np.ndarray, kappa_ref: np.ndarray
) -> tuple[AcadosModel, Function]:
    # preprocess the reference values ============================
    length = len(s_ref)
    # copy loop to beginning and end
    s_ref = np.concatenate(
        (
            -s_ref[length - 2] + s_ref[length - 81 : length - 2],
            s_ref,
            s_ref[length - 1] + s_ref[1:length],
        )
    )
    kappa_ref = np.concatenate(
        (kappa_ref[length - 80 : length - 1], kappa_ref, kappa_ref[1:length])
    )

    # compute spline interpolation for the curvature kappa =======
    kapparef_s = interpolant("kapparef_s", "bspline", [s_ref], kappa_ref)

    # set up states & controls
    s = sym_t.sym("s")
    n = sym_t.sym("n")
    alpha = sym_t.sym("alpha")
    v = sym_t.sym("v")
    T = sym_t.sym("T")
    delta = sym_t.sym("delta")
    x = vertcat(s, n, alpha, v, T, delta)

    # controls
    Tdot = sym_t.sym("Tdot")
    deltadot = sym_t.sym("deltadot")
    u = vertcat(Tdot, deltadot)

    # parameters
    p = vertcat([])

    # dynamics
    beta = C1 * delta
    F_x = (Cm1 - Cm2 * v) * T - Cr2 * v * v - Cr0 * tanh(Cr3 * v)
    a_long = F_x / m
    tpr = v * sin(beta) / l_R
    a_lat = -a_long * sin(beta) + v * tpr
    sdot_expr = (v * cos(alpha + beta)) / (1 - kapparef_s(s) * n)

    # Define model
    model = AcadosModel()
    xdot = sym_t.sym("xdot", nx)
    f_expl = vertcat(
        sdot_expr,
        v * sin(alpha + beta),
        v * C2 * delta - kapparef_s(s) * sdot_expr,
        a_long * cos(beta),
        (Tdot - T) / tau_T,
        (deltadot - delta) / tau_delta,
    )
    model.xdot = xdot
    model.f_impl_expr = xdot - f_expl
    model.f_expl_expr = f_expl

    model.x = x
    model.u = u
    model.p = p
    model.name = model_name
    model.con_h_expr = vertcat(
        a_long,
        a_lat,
        n - 0.5 * L * sin(fabs(alpha)) + 0.5 * W * cos(alpha),
        -n + 0.5 * L * sin(fabs(alpha)) + 0.5 * W * cos(alpha),
    )

    ocp = AcadosOcp()
    ocp.model = model
    ocp.dims.N = Nf
    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Qe

    Vx = np.zeros((ny, nx))
    Vx[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx = Vx

    Vu = np.zeros((ny, nu))
    Vu[-nu:] = np.eye(nu)
    ocp.cost.Vu = Vu

    Vx_e = np.zeros((ny_e, nx))
    Vx_e[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx_e = Vx_e

    ocp.cost.zl = zl
    ocp.cost.zu = zu
    ocp.cost.Zl = Zl
    ocp.cost.Zu = Zu

    # set intial references (will be overrided eithe way so the actual value doesn't matter)
    # but need them for the dimensions
    ocp.cost.yref = np.zeros(ny)
    ocp.cost.yref_e = np.zeros(ny_e)

    # setting constraints
    ocp.constraints.lbx = np.array([alpha_min, v_min, T_min, delta_min])
    ocp.constraints.ubx = np.array([alpha_max, v_max, T_max, delta_max])
    ocp.constraints.idxbx = np.array([2, 3, 4, 5])
    ocp.constraints.idxsbx = np.array([2, 3])
    ocp.constraints.lbu = np.array([Tdot_min, deltadot_min])
    ocp.constraints.ubu = np.array([Tdot_max, deltadot_max])
    if constant_controls:
        ocp.constraints.lbu *= dt
        ocp.constraints.ubu *= dt
    ocp.constraints.idxbu = np.array([0, 1])

    ocp.constraints.lh = np.array(
        [
            along_min,
            alat_min,
            -1e3,
            -1e3,
        ]
    )
    ocp.constraints.uh = np.array(
        [
            along_max,
            alat_max,
            n_max,
            n_max,
        ]
    )
    ocp.constraints.lsh = np.zeros(nsh)
    ocp.constraints.ush = np.zeros(nsh)
    ocp.constraints.idxsh = np.array(range(nsh))

    # set intial condition (same as for yref)
    ocp.constraints.x0 = np.zeros(nx)

    # set QP solver and integration
    ocp.solver_options.tf = Nf * dt
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.nlp_solver_type = "SQP_RTI"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    # ocp.solver_options.qp_solver_iter_max = 100
    ocp.solver_options.hpipm_mode = "SPEED_ABS"
    ocp.solver_options.integrator_type = "ERK"
    ocp.code_export_directory = model_name + "_gen_code"

    # create solver
    solver = AcadosOcpSolver(ocp, json_file=f"{model_name}_ocp.json", verbose=False)
    return (
        solver,
        Function("a_lat", [x, u], [a_lat]),
        Function("a_long", [x, u], [a_long]),
    )
