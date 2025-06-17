# Solvers

This repository contains tests with two function solvers: `Newton-Raphson` and `Brent`.
Both are implemented in Python using only the `math` library.


## Given equation:

All the examples here are developed around the following function.
However, the implementation are general.
Thus, you can change write your own function (and its derivative).

$$
f(R_{eq}) = -\left( \frac{t}{C \cdot \ln\left( \frac{v - I_h R_{eq}}{v_0 - I_h R_{eq}} \right)} \right) - R_{eq}
$$


## Derivative

To find the derivative of the given equation with respect to $R_{eq}$, we should carefully analyze and differentiate the expression step by step.


### Step 1: Simplify the Expression

Let's rewrite the equation for clarity:

$$
f(R_{eq}) = -\frac{t}{C} \cdot \left( \ln\left( \frac{v - I_h R_{eq}}{v_0 - I_h R_{eq}} \right) \right)^{-1} - R_{eq}
$$

Let $u = \frac{v - I_h R_{eq}}{v_0 - I_h R_{eq}}$. **Remember** $u$ is a function of $R_eq$.

$$
f(R_{eq}) = -\frac{t}{C} \cdot \left( \ln\left( u \right) \right)^{-1} - R_{eq}
$$



### Step 2: Differentiate Term by Term
We will differentiate each term with respect to $R_{eq}$.

#### Term 1: $f_1(R_{eq}) = -\frac{t}{C} \cdot \left( \ln\left( u \right) \right)^{-1}$

Using the chain rule:

$$
\frac{d}{dR_{eq}}  f_1(R_{eq}) = \frac{t}{C} \cdot \left( \ln u \right)^{-2} \cdot \frac{1}{u} \cdot \frac{du}{dR_{eq}}
$$

Now, compute $\frac{du}{dR_{eq}}$:

$$
u = \frac{v - I_h R_{eq}}{v_0 - I_h R_{eq}}
$$

$$
\frac{du}{dR_{eq}} = \frac{(-I_h)(v_0 - I_h R_{eq}) - (v - I_h R_{eq})(-I_h)}{(v_0 - I_h R_{eq})^2} = \frac{-I_h v_0 + I_h^2 R_{eq} + I_h v - I_h^2 R_{eq}}{(v_0 - I_h R_{eq})^2} = \frac{I_h (v - v_0)}{(v_0 - I_h R_{eq})^2}
$$

Substituting back:

$$
\frac{d}{dR_{eq}} f_1(R_{eq}) = \frac{t}{C} \cdot \frac{1}{(\ln u)^2} \cdot \frac{1}{u} \cdot \frac{I_h (v - v_0)}{(v_0 - I_h R_{eq})^2}
$$

Simplify $\frac{1}{u} \cdot \frac{1}{(v_0 - I_h R_{eq})^2}$:

$$
\frac{1}{u} = \frac{v_0 - I_h R_{eq}}{v - I_h R_{eq}}
$$

$$
\frac{1}{u} \cdot \frac{1}{(v_0 - I_h R_{eq})^2} = \frac{1}{(v - I_h R_{eq})(v_0 - I_h R_{eq})}
$$

Thus:

$$
\frac{d}{dR_{eq}} f_1(R_{eq}) = \frac{t I_h (v - v_0)}{C (v - I_h R_{eq})(v_0 - I_h R_{eq}) (\ln u)^2}
$$

#### Term 2: $-R_{eq}$
$$
\frac{d}{dR_{eq}} (-R_{eq}) = -1
$$

### Final Derivative:
Combine the derivatives of both terms:
$$
\frac{df}{dR_{eq}} = \frac{t I_h (v - v_0)}{C (v - I_h R_{eq})(v_0 - I_h R_{eq}) \left( \ln \left( \frac{v - I_h R_{eq}}{v_0 - I_h R_{eq}} \right) \right)^2} - 1
$$

### Simplified Form:
Let $u = \frac{v - I_h R_{eq}}{v_0 - I_h R_{eq}}$, then:
$$
\frac{df}{dR_{eq}} = \frac{t I_h (v - v_0)}{C (v - I_h R_{eq})(v_0 - I_h R_{eq}) (\ln u)^2} - 1
$$

**Important points:**
1. The term $\ln u$ appears in the denominator squared, so $u \neq 1$ (i.e., $v \neq v_0$) to avoid division by zero.
2. The derivative is valid only when $v - I_h R_{eq}$ and $v_0 - I_h R_{eq}$ are positive (since $\ln$ is defined for positive arguments).




---

# Root-Finding Techniques


# 1. Newton-Raphson

It's an iterative algorithm that uses the function and its derivative to refine an initial guess until it converges to a root.

```
def newton_raphson(f, df, x0, tol=1e-6, max_iter=100):
    """
        f: function
        df: derivative of the function
        x0: initial guess
    """
    for _ in range(max_iter):
        x1 = x0 - f(x0) / df(x0)
        if abs(x1 - x0) < tol:
            return x1
        x0 = x1
    raise ValueError("No convergence")
```

**References:**
- [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method)
- [Newton Raphson Method](https://brilliant.org/wiki/newton-raphson-method/)


# 2. Advanced Techniques

Brent's and Müller's methods are numerical algorithms for finding roots of equations, offering improvements over basic methods like bisection or Newton-Raphson.

**Comparison: Brent vs. Müller**

| Feature               | Brent's Method                          | Müller's Method                     |
|-----------------------|----------------------------------------|-------------------------------------|
| **Convergence**       | Superlinear (near quadratic)           | Quadratic (for simple roots)        |
| **Bracketing Needed?**| Yes                                    | No                                  |
| **Handles Derivatives?**| No (derivative-free)                 | No (derivative-free)               |
| **Complex Roots?**    | No                                     | Yes                                 |
| **Robustness**        | Very robust (always converges)        | May diverge with bad initial points |
| **Best For**          | Reliable, bracketed root-finding      | Polynomials/complex roots           |


### **When to Use Which Method?**
- **Use Brent's method** if:
  - You have a bracketing interval \([a, b]\).
  - You want guaranteed convergence (e.g., in scientific computing).
  - Implemented in Scipy [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brent.html).

- **Use Müller's method** if:
  - You need to find complex roots.
  - No bracketing is available (e.g., for polynomials).
  - You can find an implementation using only standard python libraries [here](https://github.com/fgittins/mullerpy).




## a. Brent's Method

Brent's method is a **hybrid algorithm** that combines **Bisection** (guaranteed convergence), **Secant method** (fast superlinear convergence), and **Inverse quadratic interpolation** (efficient when near the root).
It is robust (i.e., always converges) and efficient (i.e., fast for well-behaved functions).

For this method to work, you need to know an interval $[a, b]$ where $f(a) \cdot f(b) < 0$ (Bracketing requirement).
The method is adaptive switcing between bisection and faster methods based on behavior.
Its convergence is superlinear (close to quadratic) for smooth functions.
One advantage is that it does not need for the derivatives.


### How It Works
1. **Start with \([a, b]\)** where \(f(a)\) and \(f(b)\) have opposite signs.
2. **Try inverse quadratic interpolation** (IQI):
   - Fit a quadratic curve through 3 points and estimate the root.
   - If IQI gives a point within bounds and improves convergence, use it.
3. **Fallback to secant method** if IQI fails.
4. **Use bisection** if the above methods don't make sufficient progress.
5. **Repeat** until \(|f(x)| < \text{tolerance}\).


**Disadvantages**
- **Slightly more complex** than bisection or Newton.
- **Requires bracketing** (unlike Newton's method).



## b. Müller's Method


This method is a **generalization of the secant method** that uses **quadratic interpolation** (instead of linear).
It can find **complex roots** (unlike Brent's method) and does **not require bracketing**.


### **How It Works**
1. **Start with 3 initial points** \((x_0, x_1, x_2)\).
2. **Fit a quadratic polynomial** \(P(x)\) through these points.
3. **Find the root of \(P(x)\) closest to \(x_2\)** (using quadratic formula).
4. **Update points**: Discard the oldest point, keep the two newest.
5. **Repeat** until \(|f(x)| < \text{tolerance}\).

### **Advantages**
- **Faster than secant method** (quadratic convergence for simple roots).
- **Works without bracketing**.
- **Can find complex roots** (useful for polynomials).

### **Disadvantages**
- **May diverge** if initial points are poorly chosen.
- **More complex per iteration** than secant.
