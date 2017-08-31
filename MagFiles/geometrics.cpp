void PERC_to_normal_coordinat_transform(double *position_vector);
void PERC_to_normal_coordinat_transform(double *position_vector, double *position_vector_2);
void rotate_vector_in_3d(double *a_vector, double *rotation_center, double *rotation_axis, double alpha, double *result);
void rotate_vector_in_3d(double *a_vector, double *rotation_center, double *rotation_axis, double alpha);
void normalize_vector(double *direction);
double vector_length(double *direction);
void crossproduct (double *A, double *B, double *C);
double innerproduct (double *A, double *B);
double angle_between_vectors(double *vector_1, double *vector_2);
double clothoid_X_fun (double clothoid_A, double clothoid_L);
double clothoid_Y_fun (double clothoid_A, double clothoid_L);


void PERC_to_normal_coordinat_transform(double *position_vector) {
// transform position vectors from the PERC system(x_p,y_p,z_p) to normal coordinate system (x_p => z, y_p = -y, z_p => x)
double buffer_x_p = position_vector[1];
position_vector[1] = position_vector[3];    // x = z_p
position_vector[2] *= -1.;                  // y = -y_p
position_vector[3] = buffer_x_p;            // z = x_p
}
void PERC_to_normal_coordinat_transform(double *position_vector, double *position_vector_2) {
// transform two position vectors from the PERC system(x_p,y_p,z_p) to normal coordinate system (x_p => z, y_p = -y, z_p => x)
double buffer_x_p = position_vector[1];
position_vector[1] = position_vector[3];    // x = z_p
position_vector[2] *= -1.;                  // y = -y_p
position_vector[3] = buffer_x_p;            // z = x_p
buffer_x_p = position_vector_2[1];
position_vector_2[1] = position_vector_2[3];    // x = z_p
position_vector_2[2] *= -1.;                  // y = -y_p
position_vector_2[3] = buffer_x_p;            // z = x_p
}

void rotate_vector_in_3d(double *a_vector, double *rotation_center, double *rotation_axis, double alpha, double *result) {
// translate everything to a descard coordinat system where the centerpoint is the origin
for (int i = 1; i <= 3; i++) {a_vector[i] -= rotation_center[i];}
normalize_vector(rotation_axis);            // =n
// rotate the vector
// R(alpha)x = n(n.x) + cos(alpha) (n cross x) cross n +sin(alpha) (n cross x)
double n_cross_x[4];
crossproduct(rotation_axis,a_vector,n_cross_x);
double n_cross_x_cross_n[4];
crossproduct(n_cross_x,rotation_axis,n_cross_x_cross_n);
for (int i = 1; i<=3; i++) {
    result[i] = rotation_axis[i]*innerproduct(rotation_axis,a_vector) + cos(alpha)*n_cross_x_cross_n[i] + sin(alpha)*n_cross_x[i];
   }
// translate everything back to the origional descard coordinat system
for (int i = 1; i <= 3; i++) {result[i] += rotation_center[i];}
double alpha_2 = angle_between_vectors(result,a_vector);        // test if the rotation was good
cout << "goal angle: " << alpha << "  angle rotated: " << alpha_2 << endl;
}

void rotate_vector_in_3d(double *a_vector, double *rotation_center, double *rotation_axis, double alpha) {
// translate everything to a descard coordinat system where the centerpoint is the origin
for (int i = 1; i <= 3; i++) {a_vector[i] -= rotation_center[i];}
normalize_vector(rotation_axis);            // =n
// rotate the vector
// R(alpha)x = n(n.x) + cos(alpha) (n cross x) cross n +sin(alpha) (n cross x)
double n_cross_x[4];
crossproduct(rotation_axis,a_vector,n_cross_x);
double n_cross_x_cross_n[4];
crossproduct(n_cross_x,rotation_axis,n_cross_x_cross_n);
double result[4];
for (int i = 1; i<=3; i++) {
    result[i] = rotation_axis[i]*innerproduct(rotation_axis,a_vector) + cos(alpha)*n_cross_x_cross_n[i] + sin(alpha)*n_cross_x[i];
  }
// translate everything back to the origional descard coordinat system
for (int i = 1; i <= 3; i++) {result[i] += rotation_center[i];}
// overwrite the starting vector with the rotated one
for (int i = 1; i <= 3; i++) {a_vector[i] = result[i];}
}

void normalize_vector(double *direction) {
// normalize the direction:
double abs_direction = vector_length(direction);    // calulate the absolut value of the direction vector
for (int i = 1; i <=3 ; i++) {direction[i] /= abs_direction; }  // normalize the direction vector
}

double vector_length(double *direction) {
// calculate the length of a vector:
double abs_direction = 0.;
for (int i = 1; i <=3 ; i++) {abs_direction += pow(direction[i],2.); }
abs_direction = sqrt(abs_direction); // calulate the absolut value of the direction vector
return abs_direction;           // return the length of the vector
}

void crossproduct (double *A, double *B, double *C) {
    // AxB = C calculation
C[1] = A[2]*B[3] - A[3]*B[2];
C[2] = A[3]*B[1] - A[1]*B[3];
C[3] = A[1]*B[2] - A[2]*B[1];
}

double innerproduct (double *A, double *B) {
    // calculate the inner product of two vectors
double product = 0.;
for (int i = 1; i <= 3; i++ ) {product += A[i]*B[i];}
return product;
}

double angle_between_vectors(double *vector_1, double *vector_2) {
// cos alpha = a*b/|a|/|b|
double alpha = innerproduct(vector_1,vector_2);     // angle between vectors [rad]
alpha /= vector_length(vector_1);
alpha /= vector_length(vector_2);
alpha =  acos(alpha);
return alpha;
}


double clothoid_X_fun (double clothoid_A, double clothoid_L) {
    // calculate the x coordinate of a klothoid
  double L_A = clothoid_L/clothoid_A;
  double clothoid_X = clothoid_L + clothoid_A*( -pow(L_A,5.)/40. + pow(L_A,9.)/3456. -  pow(L_A,13.)/599040. + pow(L_A,17.)/175472640. );
  return clothoid_X;
}
double clothoid_Y_fun (double clothoid_A, double clothoid_L) {
    // calculate the x coordinate of a klothoid
  double L_A = clothoid_L/clothoid_A;
  double clothoid_Y = clothoid_A*( pow(L_A,3.)/6. - pow(L_A,7.)/336. + pow(L_A,11.)/42240. - pow(L_A,15.)/9676800. + pow(L_A,19.)/3530096640. );
  return clothoid_Y;
}
