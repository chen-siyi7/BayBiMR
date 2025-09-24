Rcpp::sourceCpp(code = '
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
Rcpp::List baymix_update_latents_cpp(
  const arma::vec& bX,
  const arma::vec& bY,
  const arma::vec& uX,
  const arma::vec& uY,
  const double theta_XY,
  const double theta_YX,
  const arma::imat& zw,        // p x 2, col0=z, col1=w (integers)
  const arma::mat& J4          // 4x4 precision of slab
) {
  const int p = bX.n_elem;

  arma::vec x(p), g(p), y(p), d(p);

  const double tXY2 = theta_XY * theta_XY;
  const double tYX2 = theta_YX * theta_YX;

  // precompute (scalar*vector uses *, NOT %)
  arma::vec q11g = uX + tXY2 * uY;
  arma::vec q12g = theta_XY * uY;
  arma::vec q22g = uY;
  arma::vec h1g  = uX % bX + theta_XY * (uY % bY);
  arma::vec h2g  = uY % bY;

  arma::vec q11d = uY + tYX2 * uX;
  arma::vec q12d = theta_YX * uX;
  arma::vec q22d = uX;
  arma::vec h1d  = uY % bY + theta_YX * (uX % bX);
  arma::vec h2d  = uX % bX;

  arma::mat Prec, JAA, R;
  arma::vec h, mu, znoise, draw;

  for (int j = 0; j < p; ++j) {
    const int zj = zw(j,0);
    const int wj = zw(j,1);

    if (zj==1 && wj==1) {
      // 4D (x, gamma, y, delta)
      Prec = J4;
      Prec(0,0) += q11g(j); Prec(0,1) += q12g(j); Prec(1,0) += q12g(j); Prec(1,1) += q22g(j);
      Prec(2,2) += q11d(j); Prec(2,3) += q12d(j); Prec(3,2) += q12d(j); Prec(3,3) += q22d(j);

      h.set_size(4);
      h(0)=h1g(j); h(1)=h2g(j); h(2)=h1d(j); h(3)=h2d(j);

      R  = arma::chol(Prec);                                    // Prec = R.t()*R
      mu = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), h));
      znoise.set_size(4); znoise.randn();
      draw = mu + arma::solve(arma::trimatu(R), znoise);

      x(j)=draw(0); g(j)=draw(1); y(j)=draw(2); d(j)=draw(3);

    } else if (zj==1 && wj==0) {
      // 3D (x, gamma, y)
      JAA = J4.submat(uvec{0,1,2}, uvec{0,1,2});
      Prec = JAA;
      Prec(0,0) += q11g(j); Prec(0,1) += q12g(j); Prec(1,0) += q12g(j); Prec(1,1) += q22g(j);
      Prec(2,2) += q11d(j);

      h.set_size(3);
      h(0)=h1g(j); h(1)=h2g(j); h(2)=h1d(j);

      R  = arma::chol(Prec);
      mu = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), h));
      znoise.set_size(3); znoise.randn();
      draw = mu + arma::solve(arma::trimatu(R), znoise);

      x(j)=draw(0); g(j)=draw(1); y(j)=draw(2); d(j)=0.0;

    } else if (zj==0 && wj==1) {
      // 3D (x, y, delta)
      JAA = J4.submat(uvec{0,2,3}, uvec{0,2,3});
      Prec = JAA;
      Prec(0,0) += q11g(j);
      Prec(1,1) += q11d(j); Prec(1,2) += q12d(j); Prec(2,1) += q12d(j); Prec(2,2) += q22d(j);

      h.set_size(3);
      h(0)=h1g(j); h(1)=h1d(j); h(2)=h2d(j);

      R  = arma::chol(Prec);
      mu = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), h));
      znoise.set_size(3); znoise.randn();
      draw = mu + arma::solve(arma::trimatu(R), znoise);

      x(j)=draw(0); g(j)=0.0; y(j)=draw(1); d(j)=draw(2);

    } else {
      // 2D (x, y)
      JAA = J4.submat(uvec{0,2}, uvec{0,2});
      Prec = JAA;
      Prec(0,0) += q11g(j);
      Prec(1,1) += q11d(j);

      h.set_size(2);
      h(0)=h1g(j); h(1)=h1d(j);

      R  = arma::chol(Prec);
      mu = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), h));
      znoise.set_size(2); znoise.randn();
      draw = mu + arma::solve(arma::trimatu(R), znoise);

      x(j)=draw(0); g(j)=0.0; y(j)=draw(1); d(j)=0.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("x")     = x,
    Rcpp::Named("gamma") = g,
    Rcpp::Named("y")     = y,
    Rcpp::Named("delta") = d
  );
}
')
