# Ground Trunth Expression for the Governing Equations

## Thomas System

### For ARGOS Families

```r
true_terms_1 <- list("sin_x2", "x1")
true_terms_2 <- list("sin_x3", "x2")
true_terms_3 <- list("sin_x1", "x3")
```

```r
true_terms_1 <- c("$sin(x_{2})$", "$x_{1}$"),
true_terms_2 <- c("$x_{2}$", "$sin(x_{3})$"),
true_terms_3 <- c("$sin(x_{1})$", "$x_{3}$"),
```

### For SINDy Families

```r
true_terms_1_py <- list("sin(1 x1)", "x0")
true_terms_2_py <- list("sin(1 x2)", "x1")
true_terms_3_py <- list("sin(1 x0)", "x2")
```

```r
    true_terms_1 <- c("$sin(x_{1})$", "$x_{0}$"),
    true_terms_2 <- c("$x_{1}$", "$sin(x_{2})$"),
    true_terms_3 <- c("$sin(x_{0})$", "$x_{2}$"),
```

## Lorenz System

### For ARGOS Families

```r
true_terms_1 <- list("x1", "x2")
true_terms_2 <- list("x1", "x1x3", "x2")
true_terms_3 <- list("x1x2", "x3")
```

```r
true_terms_1 <- c("$x_{1}$", "$x_{2}$"),
true_terms_2 <- c("$x_{1}$", "$x_{1}x_{3}$", "$x_{2}$"),
true_terms_3 <- c("$x_{1}x_{2}$", "$x_{3}$"),
```

### For SINDy Families

```r
true_terms_1_py <- list("x0", "x1")
true_terms_2_py <- list("x0", "x0 x2", "x1")
true_terms_3_py <- list("x0 x1", "x2")
```

```r
true_terms_1 <- c("$x_{0}$", "$x_{1}$"),
true_terms_2 <- c("$x_{0}$", "$x_{0}x_{2}$", "$x_{1}$"),
true_terms_3 <- c("$x_{0}x_{1}$", "$x_{2}$"),
```

## dadras System

### For ARGOS Families

```r
true_terms_1 <- list("x2", "x1", "x2x3")
true_terms_2 <- list("x2", "x1x3", "x3")
true_terms_3 <- list("x1x2", "x3")
```

```r
true_terms_1 <- c("$x_{2}$", "$x_{1}$", "$x_{2}x_{3}$"),
true_terms_2 <- c("$x_{2}$", "$x_{1}x_{3}$",  "$x_{3}$"),
true_terms_3 <- c("$x_{1}x_{2}$", "$x_{3}$"),
```

### For SINDy Families

```r
true_terms_1_py <- list("x1", "x0", "x1 x2")
true_terms_2_py <- list("x1", "x0 x2", "x2")
true_terms_3_py <- list("x0 x1", "x2")
```

```r
true_terms_1 <- c("$x_{1}$", "$x_{0}$", "$x_{1}x_{2}$"),
true_terms_2 <- c("$x_{1}$", "$x_{0}x_{2}$",  "$x_{2}$"),
true_terms_3 <- c("$x_{0}x_{1}$", "$x_{2}$"),
```

# Systems included in the main text

- lorenz
- thomas
- rossler
- maxwell-bloch
- halvorsen
- dadras []

# Systems included in the appendix

- chenlee
- linear_3d (rerun this)
