import pandas as pd
import scipy.stats as st
from sympy import *

filename = "data.csv"
p = 0.95
calculation_sample = ""

args = {}
u_a = {}
nu = {}  # degrees of freedom


def _u_a(rawdata: pd.core.series.Series, **kwargs):
    return float(rawdata.std() / sqrt(len(rawdata)))


def _import(**kwargs):
    data = pd.read_csv(filename)
    for k in data.keys():
        args.update({Symbol(k): data[k].mean()})
        u_a.update({Symbol(k): _u_a(data[k])})
        nu.update({Symbol(k): len(data[k]) - 1})
    code = ""
    for k in data.keys():
        code += k + ", "
    code = code[:-2]
    code = "global " + code + "; " + code + "= symbols(\'"
    for k in data.keys():
        code += k + " "
    exec(code + "\')")
    # print()


def y(f: Expr, **kwargs):
    return f.subs(args)


def u_y(f, append_to_sample=True, **kwargs):
    u = 0
    for x in args.keys():
        u += (diff(f, x).subs(args) * u_a[x]) ** 2  # todo

    # if append_to_sample:
    #     text = "Expansion remainder:\n"
    #     for x in args.keys():
    #         for y in args.keys():
    #             print_latex(diff(f, x, y).subs(args, hack2=True))

    return sqrt(u)


def R(f, args: dict, u, append_to_sample=True, **kwargs):
    remainder = 0
    for x in args.keys():
        for y in args.keys():
            remainder += diff(f, x, y).subs(args) * u[x] * u[y] / 2

    # if append_to_sample:
    #     text = "Expansion remainder:\n"
    #     for x in args.keys():
    #         for y in args.keys():
    #             print_latex(diff(f, x, y).subs(args, hack2=True))

    return remainder


def nu_eff(f, **kwargs):
    ret = 0
    for x in args.keys():
        ret += (diff(f, x).subs(args) * u_a[x]) ** 4 / nu[x]
    return u_y(f) ** 4 / ret


def student(nu, **kwargs):
    # linear combination of t-coefficient
    f = int(nu)
    r = nu - f
    a, b = st.t.ppf((1 + p) / 2, f), st.t.ppf((1 + p) / 2, f + 1)
    return a * (1 - r) + b * r


def U(f, **kwargs):
    if 'nu' in kwargs:
        return student(kwargs['nu']) * u_y(f, **kwargs)
    return student(nu_eff(f)) * u_y(f)


def final(f, **kwargs):
    return y(f, **kwargs), U(f, **kwargs)


def final_print(name, expr):
    name = name.replace('lambda', '_lambda')
    name = parse_expr(name, local_dict={'_lambda': Symbol('lambda')})
    print(pretty(name) + "\t = \t" + pretty(expr))


# print(R(f, args, u) , 0.8*0.0572435)
_import()
# print(U_y(m*V/t))
print(final(V))
print(final(m))
print(final(t))
print(final(V / m))
print(final(V / m, nu=2))
