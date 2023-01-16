#include <Optimizer.hpp>
#include <iostream>
#include <Const.hpp>

int main() {
    // DeformParams deformParam;
    // deformParam.B = 1.05963;
    // deformParam.C11 = 1.34514;
    // deformParam.C12 = 
    Optimizer optimizer;
    optimizer.run();
    // std::cout << "B-B: " << optimizer.getBBParams() << std::endl;
    // std::cout << "A-B: " << optimizer.getABParams() << std::endl;
    // std::cout << "A-A: " << optimizer.getAAParams() << std::endl;


    // optimizer.test();

    // Parameters params;
    // double a0 = 3.4011;
    // std::vector<Dot*> *res = new std::vector<Dot*>();
    // for (double i = 0; i < 3; ++i) {
    //     for (double j = 0; j < 3; ++j) {
    //         for (double m = 0; m < 3; ++m) {
    //             res->push_back(new Dot{i * a0, j * a0, m * a0});
    //             res->push_back(new Dot{(i + 0.5) * a0, (j + 0.5) * a0, m * a0});
    //             res->push_back(new Dot{i * a0, (j + 0.5) * a0, (m + 0.5) * a0});
    //             res->push_back(new Dot{(i + 0.5) * a0, j * a0, (m + 0.5) * a0});
    //         }
    //     }
    // }
    // params.withA0(0.086)
    // .withA1(-0.852)
    // .withP(10.96)
    // .withQ(2.278)
    // .withKsi(1.224)
    // .withCubeSide(3.615);


    //  optimizer.errorBB(params);
    
    // DeformParams deformParams = countParams(*res, eCoh, params);

    // std::cout << "Ecoh = " << eCoh  << ' '
    // << "B = " << deformParams.B << ' '
    // << "C11 = " << deformParams.C11 << ' '
    // << "C12 = " << deformParams.C12 << ' '
    // << "C44 = " << deformParams.C44 << std::endl;
    return 0;
}