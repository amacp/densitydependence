
#include <iostream>
#include <assert.h>     /* assert */
#include <Eigen/Dense>

using namespace Eigen;

MatrixXd distanceSquared(MatrixXd n1, MatrixXd n2){
    const int n1cols = n1.cols(), n1rows = n1.rows();
    const int n2cols = n2.cols(), n2rows = n2.rows();
    assert(n1cols == n2cols);
    MatrixXd temp =  (n1.array()*n1.array()).matrix()*MatrixXd::Ones(n1cols,n2rows)+ MatrixXd::Ones(n1rows,n1cols)*(n2.transpose().array()*n2.transpose().array()).matrix() - 2*n1*n2.transpose();
    return temp;
}

int main(){
    MatrixXd a(2,3); a << 1, 2, 3, 4, 5, 6;
    MatrixXd b(3,3); b << 5,6,7,3,4,5,5,4,3;
    std::cout << distanceSquared(a,b);
    return 0;
}

// + MatrixXd::Ones(n1rows,n1cols)*(n2.transpose().array()*n2.transpose().array()).matrix() - 2*n1*n2.transpose();