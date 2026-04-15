# %%
# from scipy.integrate import odeint
import re

import numpy as np


## split string 'sin(cos(x))tan(x)' to ['sin(cos(x))', 'tan(x)']
def split_funcs(string):
    return re.split(r"(?<=\))(?=[a-z])", string)


# split_funcs('sin(cos(xy))tan(x)')


def split_term_func(string):
    # string = 'sin(xz)'
    if len(re.findall("\\^", string)) == 0:
        terms = re.findall(r"(\w+)", string)
    else:
        terms = re.findall(r"(\w+\^\d+|\w+)", string)
        # terms = re.findall(r'(\w+\^?\d*|\w+\^?\d*)', string)
        # terms = re.findall(r'(\w+(?:\^\d+)?|\w+(?:\^\d+)?)', string)
        if len(terms) == 1 and len(string) > 3:
            terms1 = list(re.findall(r"(\w\d*)?(\w\d*\^\d*)|w*", terms[0])[0])
            term1_tmp = [re.findall(r"^\d\^\d", q) != [] for q in terms1]
            if any(term1_tmp) == False:
                terms = terms1
    return terms


def split_term_func2(string):
    if len(re.findall("\\^", string)) == 0:
        # Use a regular expression to capture sin, cos, etc., numbers, and variables
        terms = re.findall(r"([a-zA-Z]+|\d*\.\d+|\d+|\w)", string)
    else:
        # If there are exponents, modify the regular expression accordingly
        terms = re.findall(r"([a-zA-Z]+|\d*\.\d+|\d+|\w+\^\d+|\w+)", string)

    return terms


def split_term(string):
    # string = 'sin(xz)'
    if len(re.findall("\\^", string)) == 0:
        terms = re.findall(r"(\w\d|\w\d)", string)
        if len(terms) == 0:
            terms = re.findall(r"(\w)", string)
    else:
        terms = re.findall(r"(\w+\^\d+|\w)", string)
        # terms = re.findall(r'(.+[\^\d+]*|.+[\^\d+]*)', string)
        # terms = list(re.findall(r'(\w\d*)?(\w\d*\^\d*)', string)[0])
    return terms


# split_term('xz')
# split_term('x2^5x3^3')


## find function orders for a term e.g. find [sin, cos] from sin(cos(x))
def find_functions(string):
    # string = 'sin(cos(x))'
    # string = 'sin(cos(xy))'
    # string = 'x1^2'
    # string = 'x1^2x2'
    # string = 'x1x2^2'
    # string = 'x1^2x2^2'
    # string = 'x'
    symbols = split_term_func(string)
    basic_funs = ["sin", "cos", "tan", "log", "exp"]
    funs_order = []
    for i in basic_funs:
        try:
            index = symbols.index(i)
            funs_order.append(index)
        except ValueError:
            pass
    find_basic_funs = [symbols[i] for i in funs_order]
    var_index = [x for x in list(range(len(symbols))) if x not in funs_order]
    variables0 = [symbols[i] for i in var_index]
    if len(variables0) == 1:
        variables1 = [split_term(variables0[j]) for j in range(len(var_index))][0]
    else:
        variables1 = variables0
    # need to consider the operation order
    values1 = np.repeat(1.0, len(variables1)).tolist()
    num_index = []
    var_index = []
    for i in range(len(variables1)):
        try:
            value = eval(variables1[i])  # float(variables1[i])
            num_index.append(i)
        except:
            value = 1.0
            var_index.append(i)
        values1[i] = value
    variables = [variables1[i] for i in var_index]
    values = [values1[i] for i in num_index]
    if len(num_index) == 0:
        values = np.repeat(1.0, len(variables)).tolist()
    if len(variables) < len(values):
        variables.append("1")
    # if len(re.findall('\\^', string)) > 0:

    return (find_basic_funs, variables, values)


# find_functions('sin(cos(xy))')


def poly_order(string):
    # term = re.findall('[a-zA-Z0-9]*\\^', string)[0][:-1]
    if len(re.findall("\\^", string)) == 0:
        poly_order = 1
    else:
        poly_order = re.findall("\\^\d", string)[0][1]
    return int(poly_order)


# poly_order(split_term('x2^5x3^3')[0])


def find_which_term(string):
    return re.sub(r"\^\d+", "", string)


# find_which_term(split_term('x1^5x2^3')[0])


def basic_fun_np(string):
    if string == "sin":
        return np.sin
    if string == "cos":
        return np.cos
    if string == "tan":
        return np.tan
    if string == "log":
        return np.log
    if string == "exp":
        return np.exp


def term_comb(terms, term_names):  # input one terms, terms = [x1,x2,x3] or [x,y,z]
    # term_names = 'sin(cos(xz))y'
    # term_names = 'x1^2x2';terms=[1,2]
    # terms = [2,3,1]; term_names='x1^2'
    if term_names == "":
        return float(1)
    else:
        out_f = []
        split_terms = split_funcs(term_names)
        for k in range(len(split_terms)):
            # if len(re.findall('\\^', split_terms[k])) == 0:
            used_funcs, term_names_split, coeffs = find_functions(split_terms[k])
            poly = [
                poly_order(term_names_split[i]) for i in range(len(term_names_split))
            ]
            base_terms = [
                find_which_term(term_names_split[i])
                for i in range(len(term_names_split))
            ]
            # else:
            # poly = [int(re.search('\d',re.search('\\^\d', term_names_split[i]).group()).group()) for i in range(len(term_names_split))]

            if len(re.findall(r"\d", base_terms[0])) != 0:
                terms_index = [f"x{n+1}" for n in range(len(terms) - 1)]
                terms_index.append("t")
            else:
                if len(terms) == 3:
                    terms_index = ["x", "y", "t"]
                elif len(terms) == 4:
                    terms_index = ["x", "y", "z", "t"]
            terms_dict = {terms_index[i]: terms[i] for i in range(len(terms))}
            out = 1
            for i in range(len(term_names_split)):
                base_terms2 = base_terms[i]
                if base_terms2 in terms_dict:
                    out *= terms_dict[base_terms2] ** poly[i] * coeffs[i]
                else:
                    temp = eval(base_terms2)
                if base_terms2 == "1":  # consider 2.2*t+3.3
                    out += temp * coeffs[i]
            if len(used_funcs) > 0:
                used_funcs.reverse()
                for i in range(len(used_funcs)):
                    out = basic_fun_np(used_funcs[i])(out)
            out_f.append(out)
        return np.prod(out_f)


# term_comb([1,2,3], 'sin(cos(xz))y')

# use odeint
from scipy.integrate import odeint


def ode_eq_3d_odeint(y, t, paras, terms_name):
    n = len(y)
    y2 = np.concatenate((y, [t]))  # add t to y
    out = [[] for _ in range(n)]
    for i in range(len(terms_name)):  # i for each equation
        out[i] = 0
        for j in range(len(terms_name[i])):  # j for each term in the ith equation
            out[i] += paras[i][j] * term_comb(y2, terms_name[i][j])
    return out


def solve_ode_odeint(true_matrix, terms_name, initial, t_span):
    return odeint(ode_eq_3d_odeint, initial, t_span, args=(true_matrix, terms_name))


# new use solve_ivp
from scipy.integrate import solve_ivp


def ode_eq_3d_ivp(t, y, paras, terms_name):
    n = len(y)
    y2 = np.concatenate((y, [t]))  # add t to y
    out = np.zeros(n)
    for i in range(len(terms_name)):  # i for each equation
        for j in range(len(terms_name[i])):  # j for each term in the ith equation
            out[i] += paras[i][j] * term_comb(y2, terms_name[i][j])
    return out


def solve_ode_ivp(true_matrix, terms_name, initial, t_span, method="RK45"):
    sol = solve_ivp(
        ode_eq_3d_ivp,
        [t_span[0], t_span[-1]],
        initial,
        args=(true_matrix, terms_name),
        t_eval=t_span,
        method=method,
    )
    return sol.y.T


# %%
def generate_noisy_dynamical_systems_pyversion(
    variable_coeff: list,
    variable_names: list,
    n: int,
    dt: float,
    init_conditions: list,
    snr: int,
) -> any:
    """_summary_

    Args:
        variable_coeff (list): The coefficeints of the variables of the governing equations
        variable_names (list): The names of the variables of the governing equations
        n (int): The number of observations to generate
        dt (float): The time step of the generated observations
        init_conditions (list): Initial conditions of the dynamical system
        snr (int): snr of the generated observations

    Returns:
        any: A dataframe of the noisy dynamical system
    """
    t = np.arange(0, float(n) * dt, dt)  # np.arange(0, (n - 1) * dt, dt)
    x_t = solve_ode_odeint(variable_coeff, variable_names, init_conditions, t)
    # Convert snr (dB) to voltage
    snr_volt = 10 ** -(snr / 20)
    # Add noise (dB)
    if snr_volt != 0:
        x_init = x_t.copy()
        for i in range(int(x_t.shape[1])):
            x_t[:, i] = x_t[:, i] + snr_volt * np.random.normal(
                scale=np.std(x_init[:, i]), size=x_init[:, i].shape
            )
    return x_t


# %%
# def generate_noisy_dynamical_systems_pyversion(
#     variable_coeff: list,
#     variable_names: list,
#     n: int,
#     dt: float,
#     init_conditions: list,
#     snr: int,
# ) -> any:
#     t = np.arange(0, float(n) * dt, dt) #np.arange(0, (n - 1) * dt, dt)
#     x_t = solve_ode_odeint(variable_coeff, variable_names, init_conditions, t)
#     # Convert snr (dB) to voltage
#     snr_volt = 10 ** -(snr / 20)
#     # Add noise (dB)
#     noise_matrix = np.apply_along_axis(
#         lambda x_i: np.random.normal(0, snr_volt * np.std(x_i), len(x_i)), 1, x_t
#     )
#     xn = x_t + noise_matrix
#     return xn

# %%
# if __name__ == '__main__':
#     import matplotlib.pyplot as plt
#     from matplotlib import cm
#     ## Lorenz
#     n_obs = 2000; dt=0.01
#     t_span = np.arange(0, float(n_obs) * dt, dt)
#     true_matrix = [[10,-10],[28, -1, -1],[1, -8/3]]
#     terms_name = [['y','x'],['x','xz','y'],['xy','z']]
#     x_total = solve_ode(true_matrix,terms_name, [2,3,4], t_span)

#     plt.plot(t_span, x_total[:,0])
#     plt.plot(t_span, x_total[:,1])
#     plt.plot(t_span, x_total[:,2])
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(x_total[:,0], x_total[:,1], x_total[:,2], lw=0.5)
#     ax.set_xlabel('X Axis')
#     ax.set_ylabel('Y Axis')
#     ax.set_zlabel('Z Axis')
#     ax.set_title('Lorenz System')

#     # Lorenz parameters
#     # true_matrix_a = np.array([10, 28, -8 / 3])
#     # def lorenz_eq(x_t, t, parameters):
#     #     """Lorenz system"""
#     #     return [
#     #         parameters[0] * x_t[1] - parameters[0]*x_t[0],
#     #         parameters[1]*x_t[0]  - x_t[0] *x_t[2] - x_t[1],
#     #         x_t[0] * x_t[1] + parameters[2] * x_t[2],
#     #     ]

#     # t_span = np.arange(0, float(n_obs) * dt, dt)
#     # x_total2 = odeint(lorenz_eq, [2,3,4], t_span, args=(true_matrix_a,))
#     # plt.plot(t_span, x_total2[:,0])
#     # plt.plot(t_span, x_total2[:,1])
#     # plt.plot(t_span, x_total2[:,2])

#     ## THOMAS
#     n_obs = 2000; dt=0.1
#     t_span = np.arange(0, float(n_obs) * dt, dt)
#     true_matrix = [[1,-0.2],[1, -0.2],[1, -0.2]]
#     terms_name = [['sin(y)','x'],['sin(z)','y'],['sin(x)','z']]
#     x_total = solve_ode(true_matrix, terms_name, [1.1,1.1,-0.01], t_span)
#     # x_total = odeint(ode_eq_3d, [1.1,1.1,-0.01], t_span, args=(true_matrix_a,terms_name))

#     # plt.plot(t_span, x_total[:,0])
#     # plt.plot(t_span, x_total[:,1])
#     # plt.plot(t_span, x_total[:,2])
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(x_total[:,0], x_total[:,1], x_total[:,2], lw=0.5)
#     ax.set_xlabel('X Axis')
#     ax.set_ylabel('Y Axis')
#     ax.set_zlabel('Z Axis')
#     ax.set_title('THOMAS')

#     ## CHEN
#     n_obs = 2000; dt=0.01
#     t_span = np.arange(0, float(n_obs) * dt, dt)
#     true_matrix = [[5,-1],[-10, 1],[-0.38, 1/3]]
#     terms_name = [['x','yz'],['y','xz'],['z','xy']]
#     x_total = solve_ode(true_matrix, terms_name, [5,10,10], t_span)
#     # x_total = odeint(ode_eq_3d, [5,10,10], t_span, args=(true_matrix_a,terms_name))

#     # plt.plot(t_span, x_total[:,0])
#     # plt.plot(t_span, x_total[:,1])
#     # plt.plot(t_span, x_total[:,2])
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(x_total[:,0], x_total[:,1], x_total[:,2], lw=0.5)
#     ax.set_xlabel('X Axis')
#     ax.set_ylabel('Y Axis')
#     ax.set_zlabel('Z Axis')
#     ax.set_title('CHEN')

#     ## Lotka–Volterra
#     n_obs = 2000; dt=0.01
#     t_span = np.arange(0, float(n_obs) * dt, dt)
#     true_matrix = [[1,-1],[1, -1]]
#     terms_name = [['x','xy'],['xy','y']]
#     x_total = solve_ode(true_matrix, terms_name, [2,3], t_span)
#     # x_total = odeint(ode_eq_3d, [2,3], t_span, args=(true_matrix_a,terms_name))

#     plt.plot(x_total[:,0], x_total[:,1])

#     ## VdP
#     n_obs = 2000; dt=0.01
#     t_span = np.arange(0, float(n_obs) * dt, dt)
#     true_matrix = [[1],[1.2, -1.2, -1]]
#     terms_name = [['x2'],['x2','x1^2x2','x1']]
#     terms_name = [['y'],['y','x^2y','x']]
#     x_total = solve_ode(true_matrix, terms_name, [2,3], t_span)
#     # x_total = odeint(ode_eq_3d, [2,3], t_span, args=(true_matrix,terms_name))

#     plt.plot(x_total[:,0], x_total[:,1])


# %%
## Duffing
# n_obs = 2000; dt=0.01
# t_span = np.arange(0, float(n_obs) * dt, dt)
# true_matrix = [[1],[-1, -0.1, -5]]
# terms_name = [['x2'],['x1','x2','x1^3']]
#
# x_total = solve_ode(true_matrix, terms_name, [2,3], t_span)
# # %%
# import matplotlib.pyplot as plt
# plt.plot(x_total[:,0], x_total[:,1])

def generate_white_noise(
    n: float,
    dt: float,
    init_conditions: list,
    snr: int,
    std_scale: float = 1.0,
) -> np.ndarray:
    """Generate white noise in the same format as x_t from dynamical systems.

    Args:
        n (int): The number of observations to generate
        dt (float): The time step (for consistency, though not used in noise generation)
        init_conditions (list): Initial conditions to determine the number of variables
        snr (int): SNR of the generated noise
        std_scale (float): Standard deviation scaling factor for the noise

    Returns:
        np.ndarray: White noise array with shape (n, len(init_conditions))
    """
    # Convert n to nearest integer if it has decimal places
    n = int(round(n))

    num_variables = len(init_conditions)

    # Convert snr (dB) to voltage
    snr_volt = 10 ** -(snr / 20)

    # Generate white noise with same shape as x_t would have
    noise = np.random.normal(
        loc=0.0, scale=snr_volt * std_scale, size=(n, num_variables)
    )

    return noise
