#ifndef GRADIENT_RECONSTRUCTION_H
#define GRADIENT_RECONSTRUCTION_H

class FVMSolver;

class GradientReconstruction{
    private:
        FVMSolver* solver;
    public:
        GradientReconstruction(FVMSolver* solver);
        ~GradientReconstruction();

        // ! Reconstrução usando Least Squares...
        void reconstruct_gradients();
};

#endif