#include "phiTools.h"

constexpr float upper_limit = TMath::Pi();
constexpr float lower_limit = -TMath::Pi();

float transition_phi(float &D0del_phi)
{
    // this function reads in an angle d0del_phi and changes it so that the range is confined to (-pi/2 to 3pi/2)
    if (D0del_phi > upper_limit)
    {
        D0del_phi = D0del_phi - 2 * TMath::Pi();
    }
    if (D0del_phi < lower_limit)
    {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    if (D0del_phi < -0.5 * TMath::Pi() && D0del_phi > -1 * TMath::Pi())
    {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    return D0del_phi;
}

