#include <Dot.hpp>
#include <Const.hpp>

Dot minimizeDistance(Dot v, const Dot& period) {
    if(v.x > period.x * 2){
        v.x -= period.x;
    } else if(v.x < -period.x * 2) {
        v.x += period.x;
    }
    return v;
} 