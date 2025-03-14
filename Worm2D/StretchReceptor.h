#include "../utils.h"
class StretchReceptor {
    public:
        
        virtual void SetStretchReceptorParams(int, int, double, double) = 0;
        virtual void Update() = 0;
        virtual Params<double> getStretchReceptorParams() = 0;
};