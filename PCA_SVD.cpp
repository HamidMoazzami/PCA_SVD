/**
 * @file PCA_SVD.cpp
 * @author Hamidreza Moazzami (hamid.moazzami@gmail.com)
 * @brief  A program to compute PCA of the data matrix based on SVD and QR factorization
 * @version 0.1
 * @date 2021-12-19
 * 
 * @copyright Copyright (c) 2021 Hamidreza Moazzami
 * 
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cctype>
#include <random>
using namespace std;

/**
 * @brief A class to read the data from a file, check the validity of the values and provide a suitable input for other classes 
 * 
 */
class data_processing
{
private:
    /**
     * @brief The output matrix of the class used as the input of the other classes 
     * 
     */
    vector<vector<double>> imp_input;
    /**
     * @brief The string type of the imported data before validity checking
     * 
     */
    vector<vector<string>> wave_h_st;
    /**
     * @brief The latitude meshgrid of the grid points (used only for oceanic and atmospheric input data)
     * 
     */
    vector<vector<double>> lat_meshgrid;
    /**
     * @brief The longitude meshgrid of the grid points (used only for oceanic and atmospheric input data)
     * 
     */
    vector<vector<double>> long_meshgrid;

    /**
     * @brief Function to check the validity of the string type numbers before changing to double
     * 
     * @param str The input string type number
     * @return True when the input string is a valid number and false when is invalid
     */
    bool is_num(string str)
    {
        string str_inv;
        str_inv.resize(str.size());
        size_t n = str.size();
        for (size_t i = 0; i < n; i++)
        {
            str_inv[i] = str[n - 1];
            n--;
        }

        int64_t m = 0;
        string str_char;
        str_char.resize(1);
        for (size_t i = 0; i < str.size(); i++)
        {
            str_char = str[i];
            if (str_char.find(".") == 0)
            {
                m++;
            }
        }
        /* dot is only allowed to be in between the number and not as the first and last digit. It is checked here */
        if (str.find(".") == 0 || str_inv.find(".") == 0)
        {
            return 0;
        }
        /* Checking if there are more than one dot in the digits */
        else if (m > 1)
        {
            return 0;
        }
        else
        /* now cheching if everything else is a digit */
        {
            size_t dot_index = str.find("."); /* at this point there should be at most 1 dot in the string */
            if (dot_index != string::npos)
                str = str.replace(dot_index, 1, "0");

            bool if_all_digit = all_of(str.begin(), str.end(), ::isdigit);
            return if_all_digit;
        }
    }

public:
    /**
     * @brief The input file name
     * 
     */
    string input_file;

    /**
     * @brief Check whether the input file exist
     * 
     * @return True when exist  
     */
    bool check_file()
    {
        ifstream file(input_file);

        if (!file.is_open())
        {
            cout << "^^^Error: The"
                 << " " << input_file << " "
                 << "file cannot be opened^^^" << '\n';
            file.close();
            return 0;
        }
        file.close();

        return 1;
    }

    /**
     * @brief Reading the data from the file
     * 
     */
    void import_str()
    {

        wave_h_st.resize(0);
        ifstream file(input_file);

        string wh_st_line;
        vector<string> wh_st_line_vec;
        string wh_st;
        vector<string> wh_st_vec;

        /* Reading the data from the file line by line and save it in a the wh_st_line vector */
        while (file.good())
        {
            getline(file, wh_st_line, '\n');
            wh_st_line_vec.push_back(wh_st_line);
        }
        file.close();

        /* Creating the string matrix of the input data */
        for (size_t i = 0; i < wh_st_line_vec.size(); i++)
        {
            istringstream line_stream(wh_st_line_vec[i]);
            wh_st_vec.resize(0);

            while (!line_stream.eof())
            {
                getline(line_stream, wh_st, ',');
                wh_st_vec.push_back(wh_st);
            }
            wave_h_st.push_back(wh_st_vec);
        }
    }

    /**
     * @brief Checking the validity of created string matrix of data before changing its type to double
     * 
     * @return True when the matrix contains valid string numbers 
     */
    bool valid_num()
    {
        size_t row = wave_h_st.size();
        size_t col = wave_h_st[0].size();

        for (size_t i = 0; i < row; i++)
        {
            for (size_t j = 0; j < col; j++)
            {
                if (is_num(wave_h_st[i][j]) == 0)
                {
                    /* Printing a message when one of the values is invalid */
                    cout << "^^^Error: Invalid value in"
                         << " " << input_file << " "
                         << "file^^^ at row"
                         << " " << i + 1 << " "
                         << "and "
                         << "col"
                         << " " << j + 1 << '\n';
                    return 0;
                }
            }
        }
        return 1;
    }

    /**
     * @brief Providing a suitable matrix of type double for other classes after passing all validity checks
     * 
     */
    void import_data()
    {
        size_t row = wave_h_st.size();
        size_t col = wave_h_st[0].size();

        for (size_t i = 0; i < row; i++)
        {
            vector<double> r(col);
            imp_input.push_back(r);
        }

        for (size_t i = 0; i < row; i++)
        {
            for (size_t j = 0; j < col; j++)
            {
                imp_input[i][j] = stod(wave_h_st[i][j]);
            }
        }
    }

    /**
     * @brief Get the appropriate imported matrix of data
     * 
     * @return The matrix of data in double type
     */
    vector<vector<double>> get_imp_input()
    {
        return imp_input;
    }

    /**
     * @brief Providing the dimension of the output matrix of this class
     * 
     * @return vector of size 2 in which the first value is the row and the second is the column of the imp_input matrix 
     */
    vector<size_t> get_dim_imp_input()
    {
        vector<size_t> dim(2);
        dim[0] = imp_input.size();
        dim[1] = imp_input[0].size();
        return dim;
    }

    /**
     * @brief Finding the spatial grid of the water wave data using the latitude and longitude information. 
     * This is used for ploting of the water wave height data
     * 
     * @param latitude The latitude value of the first grid point i.e. the first element of the imp_input matrix
     * @param longitude The longitude value of the first grid point i.e. the first element of the imp_input matrix
     * @param lat_step The step size of changing the latitude in space
     * @param long_step The step size of changing the longitude in space
     */
    void lat_long(double latitude, double longitude, double lat_step, double long_step)
    {
        lat_meshgrid.resize(0);
        long_meshgrid.resize(0);

        int64_t row = get_dim_imp_input()[0];
        int64_t col = get_dim_imp_input()[1];
        vector<double> long_row(col);

        /* Finding the longitude values of all grid points */
        long_row[0] = longitude;
        for (int64_t i = 1; i < col; i++)
        {
            long_row[i] = long_row[i - 1] + long_step;
        }

        /* Creating the meshgrid of longitude values */
        for (int64_t i = 0; i < row; i++)
        {
            long_meshgrid.push_back(long_row);
        }

        /* Shaping the meshgrid of latitude values */
        vector<double> lat_row(col);
        for (int64_t i = 0; i < row; i++)
        {
            lat_meshgrid.push_back(lat_row);
        }

        /* Finding the latitude values of all the grid points and creating their meshgrid */
        double new_lat = latitude;
        for (int64_t i = 0; i < row; i++)
        {
            for (int64_t j = 0; j < col; j++)
            {
                lat_meshgrid[i][j] = new_lat;
            }
            new_lat = new_lat + lat_step;
        }
    }

    /**
     * @brief Get the meshgrid of longitude values
     * 
     * @return matrix of longitude meshgrid 
     */
    vector<vector<double>> get_long_meshgrid()
    {
        return long_meshgrid;
    }

    /**
     * @brief Get the meshgrid of latitude values
     * 
     * @return matrix of latitude meshgrid
     */
    vector<vector<double>> get_lat_meshgrid()
    {
        return lat_meshgrid;
    }

    /**
     * @brief Get the dimension of meshgrids
     * 
     * @return  vector of size 2 in which the first value is the row and the second is the column of the meshgrid matrixes  
     */
    vector<size_t> get_dim_meshgrid()
    {
        vector<size_t> grid_dim(2);
        grid_dim[0] = long_meshgrid.size();
        grid_dim[1] = long_meshgrid[0].size();
        return grid_dim;
    }
};

/**
 * @brief A class for defining a matrix and doing different matrix operations 
 * 
 */
class matrix
{
private:
    /**
     * @brief Matrix transpose result
     * 
     */
    vector<vector<double>> mat_trans;

    /**
     * @brief Matrix multiplication result
     * 
     */
    vector<vector<double>> a_times_b;

public:
    /**
     * @brief A member variable in matrix form (2D vector) to define and creat matrix whenver is needed
     * 
     */
    vector<vector<double>> mat;

    /**
     * @brief Function to shape mat matrix by providing desired dimension
     * 
     * @param row The row of mat matrix
     * @param col The column of mat matrix
     */
    void new_matrix(int64_t row, int64_t col)
    {
        vector<double> temp(col);
        for (int64_t i = 0; i < row; i++)
        {
            mat.push_back(temp);
        }
    }

    /**
     * @brief Function to transpose matrix
     * 
     * @param input_mat The input matrix to apply transposing on 
     */
    void transpose_mat(vector<vector<double>> input_mat)
    {
        mat_trans.resize(0);
        size_t row_mat_trans = input_mat[0].size();
        size_t col_mat_trans = input_mat.size();
        vector<double> vec_temp;
        vec_temp.resize(col_mat_trans);
        for (size_t i = 0; i < row_mat_trans; i++)
        {
            mat_trans.push_back(vec_temp);
        }

        for (size_t i = 0; i < input_mat.size(); i++)
        {
            for (size_t j = 0; j < input_mat[i].size(); j++)
            {
                mat_trans[j][i] = input_mat[i][j];
            }
        }
    }

    /**
     * @brief Function to do matrix multiplication
     * 
     * @param mat_a The first matrix in multiplication
     * @param mat_b The second matrix in multiplication
     */
    void matrix_mult(vector<vector<double>> mat_a, vector<vector<double>> mat_b)
    {
        a_times_b.resize(0);
        vector<double> temp_vec;
        temp_vec.resize(mat_b[0].size());
        for (size_t i = 0; i < mat_a.size(); i++)
        {
            a_times_b.push_back(temp_vec);
        }
        for (size_t i = 0; i < mat_a.size(); i++)
        {
            for (size_t k = 0; k < mat_b[0].size(); k++)
            {
                double count = 0.00;
                double mult = 0.00;
                for (size_t j = 0; j < mat_b.size(); j++)
                {
                    mult = mat_a[i][j] * mat_b[j][k];
                    count = count + mult;
                }
                a_times_b[i][k] = count;
            }
        }
    }

    /**
     * @brief Get the transposed matrix result
     * 
     * @return  The transposed matrix
     */
    vector<vector<double>> get_mat_trans()
    {
        return mat_trans;
    }

    /**
     * @brief Get the multiplication matrix result
     * 
     * @return The matrix of multiplication  
     */
    vector<vector<double>> get_mat_mult()
    {
        return a_times_b;
    }
};

/**
 * @brief It is a class to compute QR factorization. QR factorization decomposes matrix x into Q(unitary orthogonal matrix) and R(upper triangle matrix).
 It is an essential pre-requisite to find the Singular Value Decomposition (SVD) of the matrix x
 * 
 */
class qr_factor
{
protected:
    /**
     * @brief The orthogonal matrix of QR factorization
     * 
     */
    vector<vector<double>> q;
    /**
     * @brief The upper triangle matrix of QR factorization
     * 
     */
    vector<vector<double>> r;
    /**
     * @brief Row of the input matrix
     * 
     */
    size_t row;
    /**
     * @brief Column of the input matrix
     * 
     */
    size_t col;

public:
    /**
     * @brief Gram_schmidt process to compute QR factorization
     * 
     * @param input_x The input matrix which is aimed to apply QR on
     */
    void gram_schmidt(vector<vector<double>> input_x)
    {
        q.resize(0);
        r.resize(0);

        row = input_x.size();
        col = input_x[0].size();

        /* Initializing Q(row * col) matrix */
        vector<double> q_temp;
        q_temp.resize(col);
        for (size_t i = 0; i < row; i++)
        {
            q.push_back(q_temp);
        }

        /* Initializing R(col * col) matrix */
        vector<double> r_temp;
        r_temp.resize(col);
        for (size_t i = 0; i < col; i++)
        {
            r.push_back(r_temp);
        }

        /* v= Orthogonal basis vectors of the matrix which are find by projection of input_x vectors */
        /* Initialization of v */
        vector<double> v;
        v.resize(row);

        /* Finding v vectors, q and r matrixes */
        for (size_t j = 0; j < col; j++)
        {
            /* step 1 of the Gram-Schmidt process: v1 = 1st column of the input_x */
            for (size_t i = 0; i < row; i++)
            {
                v[i] = input_x[i][j];
            }
            /* Skipping the first column */
            if (j > 0)
            {
                for (size_t k = 0; k < j; k++)
                {
                    /* Find the r elements by inner product */
                    r[k][j] = 0;
                    for (size_t l = 0; l < row; l++)
                    {
                        r[k][j] = r[k][j] + q[l][k] * input_x[l][j];
                    }
                    /* Subtract the projection from v which causes v to become perpendicular to all columns of Q */
                    for (size_t l = 0; l < row; l++)
                    {
                        v[l] = v[l] - r[k][j] * q[l][k];
                    }
                }
            }
            /* Find the L2 norm of the diagonal of R */
            /* r_diag_sq is the square of R diagonal values */
            double r_diag_sq = 0;
            for (size_t l = 0; l < row; l++)
            {
                r_diag_sq = r_diag_sq + v[l] * v[l];
            }
            r[j][j] = sqrt(r_diag_sq);
            /* The orthogonal vactors are found and stored in the each column of q matrix */
            for (size_t l = 0; l < row; l++)
            {
                q[l][j] = v[l] / r[j][j];
            }
        }
    }

    /**
     * @brief Get the orthogonal matrix q of QR factorization
     * 
     * @return matrix q 
     */
    vector<vector<double>> get_q()
    {
        return q;
    }

    /**
     * @brief Get the r matrix of QR factorization
     * 
     * @return matrix r
     */
    vector<vector<double>> get_r()
    {
        return r;
    }

    /**
     * @brief Get the row size of the input matrix
     * 
     * @return row size of the input matrix 
     */
    int64_t get_row()
    {
        return row;
    }

    /**
     * @brief Get the column size of the input matrix
     * 
     * @return column of the input matrix
     */
    int64_t get_col()
    {
        return col;
    }

    /**
     * @brief Computing the dimension of the q matrix
     * 
     * @return vector of size 2 containing the dimension of q. the first element is the row and the second is the column size 
     */
    vector<size_t> get_dim_q()
    {
        vector<size_t> dim_q(2);
        dim_q[0] = q.size();
        dim_q[1] = q[0].size();
        return dim_q;
    }

    /**
     * @brief Computing the dimension of the r matrix
     * 
     * @return vector of size 2 containing the dimension of r. the first element is the row and the second is the column size 
     */
    vector<size_t> get_dim_r()
    {
        vector<size_t> dim_r(2);
        dim_r[0] = r.size();
        dim_r[1] = r[0].size();
        return dim_r;
    }
};

/**
 * @brief Computing Singular Value Decomposition (SVD). 
 Using class inheritance: qr_factor class is a base class for the svd class 
 * 
 */
class svd : public qr_factor
{
protected:
    /**
     * @brief Calling an object from the matrix class
     * 
     */
    matrix my_mat;

    /**
     * @brief Matrix of basis vectors (eigenvectors). Each column is an orthonormal basis showing one feature of the study system
    (the input matrix). The first column is the most important feature vector and the second column shows the second essential one
    and so fourth.
     * 
     */
    vector<vector<double>> u;

    /**
     * @brief A diagonal matrix in which the diagonal values show the eigenvalues correspond to the eigenvectors in u matrix. 
    They are ordered from the highest eigenvalue to lowest along the sigma matrix diameter 
     * 
     */
    vector<vector<double>> sigma;

    /**
     * @brief The third matrix in svd decomposition. v contains the eigenvectors of the tanspose(x) * x matrix 
     * 
     */
    vector<vector<double>> v;

    vector<vector<double>> x_hat;

public:
    /**
     * @brief Function to decompose matrix x into u, sigma and v matrixes
     * 
     * @param x The input matrix which is aimed to apply svd on
     */
    void decompose_svd(vector<vector<double>> x)
    {
        u.resize(0);
        sigma.resize(0);
        v.resize(0);

        /* Fiding x dimensions */
        size_t x_row = x.size();
        size_t x_col = x[0].size();

        /* Transpose of x = x_tans */
        /* Calculating x_trans by calling transpose function from matrix class */
        vector<vector<double>> x_trans;
        /* Initializing x_trans and resizing it to appropriate dimension */
        vector<double> x_trans_temp(x_row);
        for (size_t i = 0; i < x_col; i++)
        {
            x_trans.push_back(x_trans_temp);
        }
        /* Computing x_tans */
        matrix my_mat;
        my_mat.transpose_mat(x);
        x_trans = my_mat.get_mat_trans();

        /* Computing x * transpose(x):    */
        /* x_x_trans is a square matrix and its row = col = row of x */
        vector<vector<double>> x_times_x_trans;
        /* Initializing the x_x_trans matrix */
        vector<double> x_x_trans_temp(x_row);
        for (size_t i = 0; i < x_row; i++)
        {
            x_times_x_trans.push_back(x_x_trans_temp);
        }
        /*Calling matrix multiplication function from matrix class */
        my_mat.matrix_mult(x, x_trans);
        x_times_x_trans = my_mat.get_mat_mult();

        /* Computing transpose(x) * x:    */
        vector<vector<double>> x_trans_times_x;
        /* Initializing the x_trans_times_x matrix */
        vector<double> x_trans_x_temp(x_col);
        for (size_t i = 0; i < x_col; i++)
        {
            x_trans_times_x.push_back(x_trans_x_temp);
        }
        /*Calling matrix multiplication function from matrix class */
        my_mat.matrix_mult(x_trans, x);
        x_trans_times_x = my_mat.get_mat_mult();

        /* creating the matrix a and copy x * transpose(x) into it */
        vector<vector<double>> a;
        for (size_t i = 0; i < x_times_x_trans.size(); i++)
        {
            a.push_back(x_times_x_trans[i]);
        }

        /* providing appropriate size for u */
        vector<double> u_temp(x_row);
        for (size_t i = 0; i < x_row; i++)
        {
            u.push_back(u_temp);
        }

        /* Initializing u with an identity matrix */
        for (size_t i = 0; i < u.size(); i++)
        {
            u[i][i] = 1;
        }

        /* Defining a square version of sigma (x_row * x_row) and then extracting sigma (x_row * x_col) from it */
        vector<vector<double>> sigma_sq_mat;
        vector<double> sigma_temp(x_row);
        for (size_t i = 0; i < x_row; i++)
        {
            sigma_sq_mat.push_back(sigma_temp);
        }

        /* Using the inherited qr_factor functions and elements and applying it on x*transpose(x) and transpose(x)*x to find svd matrixes (u, sigma and v) */
        /* Finding u: */
        for (int64_t i = 0; i < 15; i++)
        {
            /* u and a are interatively calculated here to provide svd matrixes */
            gram_schmidt(a);
            /* By using the inherited member function (gram_schmidt), the inherited member variables q and r are computed */
            /* Using the matrix_mult function by the object from the matrix class */
            my_mat.matrix_mult(u, q);
            /* The result of the matrix multiplication is saved in "a_times_b" member variable of the matrix class 
            which is used here */
            u = my_mat.get_mat_mult();
            /* updating "a" using the member variables r and q inherited from the qr_factor class and matrix_mult member function and a_times_b 
            member variable by calling the object of the matrix class */
            my_mat.matrix_mult(r, q);
            a = my_mat.get_mat_mult();
        }

        /* Now, the updated r is a diagonal matrix. By taking the square root of diagonal values, sigma_sq_mat can be obtained */
        for (size_t i = 0; i < x_row; i++)
        {
            sigma_sq_mat[i][i] = sqrt(r[i][i]);
        }

        /* Finding sigma: */
        for (size_t i = 0; i < x_row; i++)
        {
            sigma.push_back(sigma_sq_mat[i]);
        }

        /* Fiding v: */
        /* creating the matrix b and copy transpose(x) * x into it */
        vector<vector<double>> b;
        for (size_t i = 0; i < x_trans_times_x.size(); i++)
        {
            b.push_back(x_trans_times_x[i]);
        }

        /* providing appropriate size for v */
        vector<double> v_temp(x_col);
        for (size_t i = 0; i < x_col; i++)
        {
            v.push_back(v_temp);
        }

        /* Initializing v with an identity matrix */
        for (size_t i = 0; i < v.size(); i++)
        {
            v[i][i] = 1;
        }

        /* using qr_factor and matrix class to compute v */
        /* v contains the eigenvectors of transpose(x) * x and u contains the eigenvectors of x * transpose(x) */
        for (int64_t i = 0; i < 15; i++)
        {
            /* v and "b" are interatively calculated here to provide final v */
            gram_schmidt(b);
            /* By using the inherited member function (gram_schmidt), the inherited member variables q and r are computed */
            /* Using the matrix_mult function my calling my_mat object from matrix class */
            my_mat.matrix_mult(v, q);
            /* get the multiplication result using get_mat_mult function of matrix class */
            v = my_mat.get_mat_mult();
            /* updating "b" using the member variables r and q inherited from the qr_factor class and also by calling my_mat object 
            and get_mat_mult member functions from the matrix class */
            my_mat.matrix_mult(r, q);
            b = my_mat.get_mat_mult();
        }
    }

    /**
     * @brief Compute the rank r approximation of the input matrix x
     * 
     * @param r The rank number for truncating u, sigma and v
     */
    void approximate_x(int64_t r)
    {
        x_hat.resize(0);
        vector<vector<double>> reduced_u;
        vector<vector<double>> reduced_v_trans;
        vector<vector<double>> reduced_sigma;

        /* Computing rank r of u */
        for (size_t i = 0; i < u.size(); i++)
        {
            vector<double> temp_r(r);
            reduced_u.push_back(temp_r);
        }

        for (size_t i = 0; i < u.size(); i++)
        {
            for (int64_t j = 0; j < r; j++)
            {
                reduced_u[i][j] = u[i][j];
            }
        }

        /* Computing rank r of sigma */
        for (int64_t i = 0; i < r; i++)
        {
            vector<double> temp_r2(r);
            reduced_sigma.push_back(temp_r2);
        }

        for (int64_t i = 0; i < r; i++)
        {
            reduced_sigma[i][i] = sigma[i][i];
        }

        /* Computing rank r of transpose(v) */
        my_mat.transpose_mat(v);
        for (int64_t i = 0; i < r; i++)
        {
            reduced_v_trans.push_back(my_mat.get_mat_trans()[i]);
        }

        /* Computing the approximation of x = reduced_u * reduced_sigma * reduced_v_trans */
        my_mat.matrix_mult(reduced_u, reduced_sigma);
        vector<vector<double>> mult1;
        for (size_t i = 0; i < reduced_u.size(); i++)
        {
            mult1.push_back(my_mat.get_mat_mult()[i]);
        }

        my_mat.matrix_mult(mult1, reduced_v_trans);

        for (size_t i = 0; i < u.size(); i++)
        {
            x_hat.push_back(my_mat.get_mat_mult()[i]);
        }
    }

    /**
     * @brief Get the u matrix
     * 
     * @return u matrix 
     */
    vector<vector<double>> get_u_svd()
    {
        return u;
    }

    /**
     * @brief Get the sigma matrix
     * 
     * @return sigma matrix
     */
    vector<vector<double>> get_sigma_svd()
    {
        return sigma;
    }

    /**
     * @brief Get the v matrix
     * 
     * @return v matrix 
     */
    vector<vector<double>> get_v_svd()
    {
        return v;
    }

    /**
     * @brief Get the approximated version of X
     * 
     * @return x_hat matrix 
     */
    vector<vector<double>> get_x_hat()
    {
        return x_hat;
    }

    /**
     * @brief Comput the dimension of u
     * 
     * @return vector of size 2 containing the dimension of u. The first element is the size of row and 
     the second is that of column
     */
    vector<int64_t> get_dimension_u()
    {
        vector<int64_t> dim_u(2);
        dim_u[0] = u.size();
        dim_u[1] = u[0].size();
        return dim_u;
    }

    /**
     * @brief Compute the dimension of sigma
     * 
     * @return vector of size 2 containing the dimension of sigma. The first element is the size of row and 
     the second is that of column 
     */
    vector<int64_t> get_dimension_sigma()
    {
        vector<int64_t> dim_sigma(2);
        dim_sigma[0] = sigma.size();
        dim_sigma[1] = sigma[0].size();
        return dim_sigma;
    }

    /**
     * @brief Compute the dimension of v
     * 
     * @return vector of size 2 containing the dimension of v. The first element is the size of row and 
     the second is that of column 
     */
    vector<int64_t> get_dimension_v()
    {
        vector<int64_t> dim_v(2);
        dim_v[0] = v.size();
        dim_v[1] = v[0].size();
        return dim_v;
    }
};

/**
 * @brief Computing Principal Component Analysis (PCA)
 * 
 */
class pca : public svd
{
private:
    /**
     * @brief Matrix of scores or principal components
     * 
     */
    vector<vector<double>> score;
    /**
     * @brief eigenvector matrix of pca 
     * 
     */
    vector<vector<double>> u_pca;

    /**
     * @brief Diagonal matrix of eigenvalues. each eigenvalue corresponds to a column in u_pca
     * 
     */
    vector<vector<double>> sigma_pca;

    /**
     * @brief orthonormal matrix of loadings
     * 
     */
    vector<vector<double>> loading;

    /**
     * @brief centralizing the feature vectors of the input matrix
     * 
     * @param x The input matrix which is aimed to apply pca on. Here is the matrix output of data_processing class
     * @return matrix of centeralized version of x 
     */
    vector<vector<double>> centralize_x(vector<vector<double>> x)
    {
        size_t x_row = x.size();
        size_t x_col = x[0].size();

        /* Finding the mean of each column and store them in a vector */
        vector<double> x_mean(x_col);
        for (size_t j = 0; j < x_col; j++)
        {
            double sum = 0;
            int count = 0;
            for (size_t i = 0; i < x_row; i++)
            {
                sum = sum + x[i][j];
                count++;
            }
            x_mean[j] = sum / count;
        }

        /* Subtracting the mean values from each coresponding column */
        /* x_tilde = matrix of centralized values. First, initializing it */
        vector<vector<double>> x_tilde;
        for (size_t i = 0; i < x_row; i++)
        {
            vector<double> temp(x_col);
            x_tilde.push_back(temp);
        }

        /* Finding x_tilde */
        for (size_t i = 0; i < x_row; i++)
        {
            for (size_t j = 0; j < x_col; j++)
            {
                x_tilde[i][j] = x[i][j] - x_mean[j];
            }
        }
        return x_tilde;
    }

public:
    /**
     * @brief Function to decompose x into u_pca, sigma_pca and v_pca
     * 
     * @param x Input matrix
     */
    void decompose_pca(vector<vector<double>> x)
    {
        size_t x_row = x.size();
        size_t x_col = x[0].size();

        /* Fiding centralized x result and save it */
        vector<vector<double>> x_tilde_res;
        for (size_t i = 0; i < x_row; i++)
        {
            vector<double> temp(x_col);
            x_tilde_res.push_back(temp);
        }
        x_tilde_res = centralize_x(x);

        size_t x_tilde_row = x_tilde_res.size();
        size_t x_tilde_col = x_tilde_res[0].size();

        /* Shaping and initializing u_pca */
        for (size_t i = 0; i < x_tilde_row; i++)
        {
            vector<double> temp(x_tilde_row);
            u_pca.push_back(temp);
        }

        /* Shaping and initializing sigma_pca */
        for (size_t i = 0; i < x_tilde_row; i++)
        {
            vector<double> temp(x_tilde_row);
            sigma_pca.push_back(temp);
        }

        /* Shaping and initializing loading */
        for (size_t i = 0; i < x_tilde_col; i++)
        {
            vector<double> temp(x_tilde_col);
            loading.push_back(temp);
        }

        /* Applying the inherited member function from the svd class on x_tilde_res */
        decompose_svd(x_tilde_res);
        /* Finding the pca matrixes using inherited functions from the svd class */
        u_pca = get_u_svd();
        sigma_pca = get_sigma_svd();
        loading = get_v_svd();

        /* Computing the score matrix. score = u_pca * sigma_pca */
        score.resize(0);
        my_mat.matrix_mult(u_pca, sigma_pca);

        size_t score_row = my_mat.get_mat_mult().size();

        for (size_t i = 0; i < score_row; i++)
        {
            score.push_back(my_mat.get_mat_mult()[i]);
        }
    }

    /**
     * @brief Get the matrix of eigenvectors
     * 
     * @return u_pca matrix
     */
    vector<vector<double>> get_u_pca()
    {
        return u_pca;
    }

    /**
     * @brief Get the matrix of eigenvalues 
     * 
     * @return sigma_pca matrix 
     */
    vector<vector<double>> get_sigma_pca()
    {
        return sigma_pca;
    }

    /**
     * @brief Get the orthonormal loading matrix 
     * 
     * @return loading matrix 
     */
    vector<vector<double>> get_loading()
    {
        return loading;
    }

    /**
     * @brief Get the score (matrix of principal components)
     * 
     * @return score matrix 
     */
    vector<vector<double>> get_score()
    {
        return score;
    }

    /**
     * @brief Get the dimension of the u_pca matrix
     * 
     * @return vector of size 2 containing the dimension of u_pca. The first element is the size of row and 
     the second is that of column   
     */
    vector<int64_t> get_dim_u_pca()
    {
        vector<int64_t> dim_u_pca(2);
        dim_u_pca[0] = u_pca.size();
        dim_u_pca[1] = u_pca[0].size();
        return dim_u_pca;
    }

    /**
     * @brief Get the dimension of the sigma_pca
     * 
     * @return vector of size 2 containing the dimension of sigma_pca. The first element is the size of row and 
     the second is that of column  
     */
    vector<int64_t> get_dim_sigma_pca()
    {
        vector<int64_t> dim_sigma_pca(2);
        dim_sigma_pca[0] = sigma_pca.size();
        dim_sigma_pca[1] = sigma_pca[0].size();
        return dim_sigma_pca;
    }

    /**
     * @brief Get the dimension of the loading
     * 
     * @return vector of size 2 containing the dimension of loading. The first element is the size of row and 
     the second is that of column  
     */
    vector<size_t> get_dim_loading()
    {
        vector<size_t> dim_loading(2);
        dim_loading[0] = loading.size();
        dim_loading[1] = loading[0].size();
        return dim_loading;
    }

    /**
     * @brief Get the dimension of the score
     * 
     * @return vector of size 2 containing the dimension of the score matrix. The first element is the size of row and 
     the second is that of column 
     */
    vector<size_t> get_dim_score()
    {
        vector<size_t> dim_score(2);
        dim_score[0] = score.size();
        dim_score[1] = score[0].size();
        return dim_score;
    }
};

int main()
{

    /* Importing data and validity checking */

    data_processing d1;
    d1.input_file = "Gulf_input.csv";

    /* Checking whether the file exist */
    if (!d1.check_file())
    {
        return -1;
    }

    d1.import_str();
    /* Checking the validity of the data */
    if (!d1.valid_num())
    {
        return -1;
    }

    /* getting the validated data */
    d1.import_data();

    //Applying QR factorization on the input data
    qr_factor my_qr1;
    my_qr1.gram_schmidt(d1.get_imp_input());

    /* Creating the r matrix by calling an object of the matrix class */
    matrix m1;
    m1.new_matrix(my_qr1.get_dim_r()[0], my_qr1.get_dim_r()[1]);
    m1.mat = my_qr1.get_r();

    /* Creating the q matrix */
    matrix m2;
    m2.new_matrix(my_qr1.get_dim_q()[0], my_qr1.get_dim_q()[1]);
    m2.mat = my_qr1.get_q();

    // Applying SVD on the input data
    svd my_svd1;
    my_svd1.decompose_svd(d1.get_imp_input());

    /* Creating the u matrix */
    matrix m3;
    m3.new_matrix(my_svd1.get_dimension_u()[0], my_svd1.get_dimension_u()[1]);
    m3.mat = my_svd1.get_u_svd();

    /* Creating the sigma matrix */
    matrix m4;
    m4.new_matrix(my_svd1.get_dimension_sigma()[0], my_svd1.get_dimension_sigma()[1]);
    m4.mat = my_svd1.get_sigma_svd();

    /* Creating the v matrix */
    matrix m5;
    m5.new_matrix(my_svd1.get_dimension_v()[0], my_svd1.get_dimension_v()[1]);
    m5.mat = my_svd1.get_v_svd();

    /* computing the rank 5 approximation of the input matrix */
    my_svd1.approximate_x(5);
    matrix x_hat;
    x_hat.new_matrix(my_svd1.get_x_hat().size(), my_svd1.get_x_hat()[0].size());
    x_hat.mat = my_svd1.get_x_hat();

    // Applying PCA
    pca my_pca1;
    my_pca1.decompose_pca(d1.get_imp_input());

    /* Creating the u_pca matrix */
    matrix m6;
    m6.new_matrix(my_pca1.get_dim_u_pca()[0], my_pca1.get_dim_u_pca()[1]);
    m6.mat = my_pca1.get_u_pca();

    /* Creating the sigma_pca matrix */
    matrix m7;
    m7.new_matrix(my_pca1.get_dim_sigma_pca()[0], my_pca1.get_dim_sigma_pca()[1]);
    m7.mat = my_pca1.get_sigma_pca();

    /* Creating the loading matrix */
    matrix m8;
    m8.new_matrix(my_pca1.get_dim_loading()[0], my_pca1.get_dim_loading()[1]);
    m8.mat = my_pca1.get_loading();

    /* Creating the score matrix */
    matrix score;
    score.new_matrix(my_pca1.get_dim_score()[0], my_pca1.get_dim_score()[1]);
    score.mat = my_pca1.get_score();

    // Creating the meshgrid of the water wave height input data
    d1.lat_long(29.04, 48.96, -0.04, 0.04);

    /* Creating the meshgrid matrix of longitude values by calling an object of the matrix class */
    matrix m9;
    m9.new_matrix(d1.get_dim_meshgrid()[0], d1.get_dim_meshgrid()[1]);
    m9.mat = d1.get_long_meshgrid();

    /* Creating the meshgrid matrix of latitude values by calling an object of the matrix class */
    matrix m10;
    m10.new_matrix(d1.get_dim_meshgrid()[0], d1.get_dim_meshgrid()[1]);
    m10.mat = d1.get_lat_meshgrid();

    //Printing the Outputs
    string gulf_r = "Gulf_qr_r.txt";
    string gulf_q = "Gulf_qr_q.txt";
    string gulf_u_svd = "Gulf_u_svd.txt";
    string gulf_sigma_svd = "Gulf_sigma_svd.txt";
    string gulf_v_svd = "Gulf_v_svd.txt";
    string gulf_u_pca = "Gulf_u_pca.txt";
    string gulf_sigma_pca = "Gulf_sigma_pca.txt";
    string gulf_loading = "Gulf_loading.txt";
    string gulf_score = "Gulf_score.txt";
    string gulf_latlong_meshgrid = "Gulf_latlong_meshgrid.txt";
    string gulf_x_hat = "Gulf_x_hat.txt";

    ofstream fileOut_1(gulf_r);
    if (!fileOut_1.is_open())
    {
        cout << "Error in opening the " << gulf_r << '\n';
        return -1;
    }
    for (size_t i = 0; i < m1.mat.size(); i++)
    {
        for (size_t j = 0; j < m1.mat[0].size(); j++)
        {
            fileOut_1 << " " << m1.mat[i][j];
        }
        fileOut_1 << '\n';
    }
    fileOut_1 << '\n';

    ofstream fileOut_2(gulf_q);
    if (!fileOut_2.is_open())
    {
        cout << "Error in opening the " << gulf_q << '\n';
        return -1;
    }
    for (size_t i = 0; i < m2.mat.size(); i++)
    {
        for (size_t j = 0; j < m2.mat[0].size(); j++)
        {
            fileOut_2 << " " << m2.mat[i][j];
        }
        fileOut_2 << '\n';
    }
    fileOut_2 << '\n';

    ofstream fileOut_3(gulf_u_svd);
    if (!fileOut_3.is_open())
    {
        cout << "Error in opening the " << gulf_u_svd << '\n';
        return -1;
    }
    for (size_t i = 0; i < m3.mat.size(); i++)
    {
        for (size_t j = 0; j < m3.mat[0].size(); j++)
        {
            fileOut_3 << " " << m3.mat[i][j];
        }
        fileOut_3 << '\n';
    }
    fileOut_3 << '\n';

    ofstream fileOut_4(gulf_sigma_svd);
    if (!fileOut_4.is_open())
    {
        cout << "Error in opening the " << gulf_sigma_svd << '\n';
        return -1;
    }
    for (size_t i = 0; i < m4.mat.size(); i++)
    {
        for (size_t j = 0; j < m4.mat[0].size(); j++)
        {
            fileOut_4 << " " << m4.mat[i][j];
        }
        fileOut_4 << '\n';
    }
    fileOut_4 << '\n';

    ofstream fileOut_5(gulf_v_svd);
    if (!fileOut_5.is_open())
    {
        cout << "Error in opening the " << gulf_v_svd << '\n';
        return -1;
    }
    for (size_t i = 0; i < m5.mat.size(); i++)
    {
        for (size_t j = 0; j < m5.mat[0].size(); j++)
        {
            fileOut_5 << " " << m5.mat[i][j];
        }
        fileOut_5 << '\n';
    }
    fileOut_5 << '\n';

    ofstream fileOut_5a(gulf_x_hat);
    if (!fileOut_5a.is_open())
    {
        cout << "Error in opening the " << gulf_x_hat << '\n';
        return -1;
    }
    for (size_t i = 0; i < x_hat.mat.size(); i++)
    {
        for (size_t j = 0; j < x_hat.mat[0].size(); j++)
        {
            fileOut_5a << " " << x_hat.mat[i][j];
        }
        fileOut_5a << '\n';
    }
    fileOut_5a << '\n';

    ofstream fileOut_6(gulf_u_pca);
    if (!fileOut_6.is_open())
    {
        cout << "Error in opening the " << gulf_u_pca << '\n';
        return -1;
    }
    for (size_t i = 0; i < m6.mat.size(); i++)
    {
        for (size_t j = 0; j < m6.mat[0].size(); j++)
        {
            fileOut_6 << " " << m6.mat[i][j];
        }
        fileOut_6 << '\n';
    }
    fileOut_6 << '\n';

    ofstream fileOut_7(gulf_sigma_pca);
    if (!fileOut_7.is_open())
    {
        cout << "Error in opening the " << gulf_sigma_pca << '\n';
        return -1;
    }
    for (size_t i = 0; i < m7.mat.size(); i++)
    {
        for (size_t j = 0; j < m7.mat[0].size(); j++)
        {
            fileOut_7 << " " << m7.mat[i][j];
        }
        fileOut_7 << '\n';
    }
    fileOut_7 << '\n';

    ofstream fileOut_8(gulf_loading);
    if (!fileOut_8.is_open())
    {
        cout << "Error in opening the " << gulf_loading << '\n';
        return -1;
    }
    for (size_t i = 0; i < m8.mat.size(); i++)
    {
        for (size_t j = 0; j < m8.mat[0].size(); j++)
        {
            fileOut_8 << " " << m8.mat[i][j];
        }
        fileOut_8 << '\n';
    }
    fileOut_8 << '\n';

    ofstream fileOut_8a(gulf_score);
    if (!fileOut_8a.is_open())
    {
        cout << "Error in opening the " << gulf_score << '\n';
        return -1;
    }
    for (size_t i = 0; i < score.mat.size(); i++)
    {
        for (size_t j = 0; j < score.mat[0].size(); j++)
        {
            fileOut_8a << " " << score.mat[i][j];
        }
        fileOut_8a << '\n';
    }
    fileOut_8a << '\n';

    ofstream fileOut_9(gulf_latlong_meshgrid);
    if (!fileOut_9.is_open())
    {
        cout << "Error in opening the " << gulf_latlong_meshgrid << '\n';
        return -1;
    }

    fileOut_9 << "Longitude_Meshgrid" << '\n';
    for (size_t i = 0; i < m9.mat.size(); i++)
    {
        for (size_t j = 0; j < m9.mat[0].size(); j++)
        {
            fileOut_9 << " " << m9.mat[i][j];
        }
        fileOut_9 << '\n';
    }
    fileOut_9 << '\n';

    fileOut_9 << "Latitude_Meshgrid" << '\n';
    for (size_t i = 0; i < m10.mat.size(); i++)
    {
        for (size_t j = 0; j < m10.mat[0].size(); j++)
        {
            fileOut_9 << " " << m10.mat[i][j];
        }
        fileOut_9 << '\n';
    }
    fileOut_9 << '\n';

    cout << "INFO:" << '\n';
    cout << "The output results are now ready in the following files:" << '\n';
    cout << "            " << gulf_r << ": the upper triangle matrix of QR factorization" << '\n';
    cout << "            " << gulf_q << ": the unitary matrix of QR factorization" << '\n';
    cout << "           " << gulf_u_svd << ": the matrix of eigenvectors in SVD" << '\n';
    cout << "       " << gulf_sigma_svd << ": the diagonal matrix of eigenvalues in SVD" << '\n';
    cout << "           " << gulf_v_svd << ": the third matrix in SVD which is an orthonormal matrix" << '\n';
    cout << "           " << gulf_x_hat << ": the approximated rank r of x input matrix by SVD" << '\n';
    cout << "           " << gulf_u_pca << ": the eigenvector matrix in PCA" << '\n';
    cout << "       " << gulf_sigma_pca << ": the matrix of eigenvalues in PCA" << '\n';
    cout << "         " << gulf_loading << ": the orthonormal matrix of loadings" << '\n';
    cout << "           " << gulf_score << ": the orthogonal matrix of principal components" << '\n';
    cout << gulf_latlong_meshgrid << ": the meshgrid of latitude and longitude values of input water wave height data" << '\n';
}
