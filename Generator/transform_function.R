transform_function = function(x, eps=0, b=1) {
  (1+eps) * (1+exp(b)) / (1 + exp(b-x)) - eps
}
