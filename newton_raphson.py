import math


def newton_raphson(equation, derivative, R_eq_guess, args, tol=1e-6, max_iter=100):
    # 
    # Ref.: https://en.wikipedia.org/wiki/Newton%27s_method
    R_eq = R_eq_guess
    for _ in range(max_iter):
        f_value = equation(R_eq, *args)
        f_derivative = derivative(R_eq, *args)

        if abs(f_derivative) < tol:  # Prevent division by a very small number
            break

        R_eq_new = R_eq - f_value / f_derivative

        if abs(R_eq_new - R_eq) < tol:
            return R_eq_new

        R_eq = R_eq_new

    return None  # If no solution is found