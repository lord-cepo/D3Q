#include "quadrature_general.hpp"

class Phonons
{
public:
    std::array<double, 8> values;
    std::array<std::array<double, 3>, 8> gradients;
    Phonons(double *values, double *gradients)
    {
        for (int i = 0; i < 8; i++)
        {
            this->values[i] = values[i];
            for (int j = 0; j < 3; j++)
            {
                this->gradients[i][j] = gradients[3 * i + j];
            }
        }
    }

    template <typename T>
    T operator()(const algoim::uvector<T, 3> &x) const
    {
        T res = T(0.0);
        int m = 0;
        std::array<std::array<T, 3>, 2> h;
        std::array<std::array<T, 3>, 2> h1;
        for (int i = 0; i < 3; i++)
        {
            // TODO: da rivedere
            h[0][i] = 1 - 3 * x(i) * x(i) + 2 * x(i) * x(i) * x(i);
            h[1][i] = x(i) * x(i) * (3 - 2 * x(i));
            h1[0][i] = x(i) * (1 - x(i)) * (1 - x(i));
            h1[1][i] = x(i) * x(i) * (x(i) - 1);
        }

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                {
                    m = i * 4 + j * 2 + k;
                    res += h[i][0] * h[j][1] * h[k][2] * values[m] +
                           gradients[m][0] * h1[i][0] * h[j][1] * h[k][2] +
                           gradients[m][1] * h[i][0] * h1[j][1] * h[k][2] +
                           gradients[m][2] * h[i][0] * h[j][1] * h1[k][2];
                }
        return res;
    }

    template <typename T>
    algoim::uvector<T, 3> grad(const algoim::uvector<T, 3> &x) const
    {
        algoim::uvector<T, 3> res;
        for (int i = 0; i < 3; i++)
            res(i) = T(0.0);
        int m = 0;
        std::array<std::array<T, 3>, 2> h;
        std::array<std::array<T, 3>, 2> h1;
        std::array<std::array<T, 3>, 2> h_der;
        std::array<std::array<T, 3>, 2> h1_der;
        for (int i = 0; i < 3; i++)
        {
            h[0][i] = 1 - 3 * x(i) * x(i) + 2 * x(i) * x(i) * x(i);
            h1[0][i] = x(i) * (1 - x(i)) * (1 - x(i));
            h[1][i] = x(i) * x(i) * (3 - 2 * x(i));
            h1[1][i] = x(i) * x(i) * (x(i) - 1);

            h_der[0][i] = -6 * x(i) + 6 * x(i) * x(i);
            h1_der[0][i] = 3 * x(i) * x(i) - 4 * x(i) + 1;
            h_der[1][i] = -h_der[0][i];
            h1_der[1][i] = 3 * x(i) * x(i) - 2 * x(i);
        }

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                {
                    m = i * 4 + j * 2 + k;
                    res(0) += h_der[i][0] * h[j][1] * h[k][2] * values[m] +
                              gradients[m][0] * h1_der[i][0] * h[j][1] * h[k][2] +
                              gradients[m][1] * h_der[i][0] * h1[j][1] * h[k][2] +
                              gradients[m][2] * h_der[i][0] * h[j][1] * h1[k][2];

                    res(1) += h[i][0] * h_der[j][1] * h[k][2] * values[m] +
                              gradients[m][0] * h1[i][0] * h_der[j][1] * h[k][2] +
                              gradients[m][1] * h[i][0] * h1_der[j][1] * h[k][2] +
                              gradients[m][2] * h[i][0] * h_der[j][1] * h1[k][2];

                    res(2) += h[i][0] * h[j][1] * h_der[k][2] * values[m] +
                              gradients[m][0] * h1[i][0] * h[j][1] * h_der[k][2] +
                              gradients[m][1] * h[i][0] * h1[j][1] * h_der[k][2] +
                              gradients[m][2] * h[i][0] * h[j][1] * h1_der[k][2];
                }
        return res;
    }
};

extern "C"
{
    void delta_nodes(double *values, double *gradients, int quality, int qorder,
                     int &number_of_nodes, double *&weights, double *&coords)
    {
        Phonons ph = Phonons(values, gradients);
        algoim::QuadratureRule<3> q;
        int m = 0;
        double sum = 0.0;
        std::vector<algoim::QuadratureRule<3>::Node> nodes;
        number_of_nodes = 0;
        double dx = 1.0 / quality;
        for (int i = 0; i < quality; ++i)
            for (int j = 0; j < quality; ++j)
                for (int k = 0; k < quality; ++k)
                {
                    algoim::uvector<double, 3> xmin{0.0 + i * dx, 0.0 + j * dx, 0.0 + k * dx};
                    algoim::uvector<double, 3> xmax{0.0 + i * dx + dx, 0.0 + j * dx + dx, 0.0 + k * dx + dx};
                    q = algoim::quadGen<3>(ph, algoim::HyperRectangle<double, 3>(xmin, xmax), 3, -1, qorder);
                    // std::cout << "Number of nodes: " << q.nodes.size() << std::endl;
                    m = i * quality * quality + j * quality + k;
                    nodes.insert(nodes.end(), q.nodes.begin(), q.nodes.end());
                }
        number_of_nodes = nodes.size();
        weights = new double[number_of_nodes];
        coords = new double[3 * number_of_nodes];
        for (int i = 0; i < number_of_nodes; i++)
        {
            weights[i] = nodes[i].w;
            sum += weights[i];
            for (int j = 0; j < 3; j++)
            {
                coords[3 * i + j] = nodes[i].x(j);
            }
        }
        // std::cout << "Sum of weights: " << sum << std::endl;
    }

    void free_c_ptr(double *&myPointer)
    {
        delete myPointer; // freed memory
        myPointer = NULL; // pointed dangling ptr to NULL
    }
} // extern "C"