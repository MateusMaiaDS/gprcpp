#include <RcppEigen.h>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Type defining as standard template
typedef Map<MatrixXd> MapMatd; // Mapping a double matrix
typedef Map<VectorXd> MapVecd; // Mapping a double vector
