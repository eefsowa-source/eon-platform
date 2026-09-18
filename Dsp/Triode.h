#pragma once
#include <cmath>
#include <algorithm>

namespace eon {

// ---------------------------------------------------------------------------
// Koren triode model + REAL Newton-Raphson plate-voltage solve.
//
// Koren model (12AX7-flavoured defaults):
//   E1 = (Vp/Kp) * ln(1 + exp(Kp*(1/mu + Vg/sqrt(Kvb + Vp^2))))
//   Ip = E1^ex / Kg                                        (E1 > 0)
//
// Circuit: plate supply B+ through load Rload.
//   f(Vp) = (B+ - Vp)/Rload - Ip(Vg, Vp) = 0
//
// Newton step uses the ANALYTIC dIp/dVp (the original used a fabricated
// derivative `Ip*0.01` and solved an equation that wasn't the circuit).
// Solution is warm-started from the previous sample — converges in ~2 iters.
// ---------------------------------------------------------------------------

struct TriodeStage
{
    // Koren parameters
    float mu  = 100.f;   // amplification factor
    float ex  = 1.4f;    // exponent
    float kg  = 1060.f;  // (roughly) inverse transconductance scale
    float kp  = 600.f;
    float kvb = 300.f;

    // Circuit
    float Bplus  = 250.f;     // HT supply [V]
    float Rload  = 100000.f;  // total plate+cathode load [ohm]
    float biasVg = -1.5f;     // DC grid bias [V]
    float inVolts = 3.0f;     // input sensitivity: full-scale -> volts of Vg swing

    int iterations = 2;

    // --- optional self-biasing cathode (Rk || Ck to ground) ---------------
    // Off by default: cathodeRk == 0 keeps cathodeVk == 0 and the output is
    // bit-identical to the grounded-cathode version.
    //
    // Enabled via setCathode(): the cathode rides on Rk||Ck so the average
    // plate current lifts Vk, which reduces Vgk and self-biases the tube
    // colder — program-dependent "cathode sag" compression (ported from the
    // Tube comp TriodeStage, which uses the same mechanism for vari-mu feel).
    float  cathodeRk   = 0.f;     // cathode resistor [ohm], 0 = grounded
    float  cathodeCk   = 25e-6f;  // bypass capacitor [F]
    double cathodePole = 0.0;     // exp(-1/(Rk*Ck*rate)) — set by setCathode
    double cathodeVk   = 0.0;     // tracked cathode voltage [V]
    float  lastIp      = 0.f;     // plate current from last solve [A]

    // state
    float VpPrev = 150.f;
    float Vq = 150.f;         // quiescent plate voltage (solved at reset)
    float outScale = 1.f;

    // Fit the Rk||Ck cathode network. rateHz is the *processing* rate —
    // under oversampling pass the oversampled rate, not the host rate.
    // Call before reset() so the quiescent solve sees the fitted network.
    void setCathode (float rkOhms, float ckFarads, double rateHz)
    {
        cathodeRk = rkOhms;
        cathodeCk = ckFarads;
        cathodePole = rkOhms > 0.f
            ? std::exp (-1.0 / ((double) rkOhms * ckFarads * rateHz)) : 0.0;
    }

    void reset()
    {
        VpPrev = 150.f;
        cathodeVk = 0.0;
        lastIp = 0.f;
        // Quiescent point. With a self-bias cathode the operating point is
        // itself signal-dependent: iterate plate-solve <-> Vk = Ip*Rk until
        // they agree (same fixed-point seeding as the Tube comp stage).
        const int n = (cathodeRk > 0.f && cathodePole > 0.0) ? 8 : 1;
        for (int i = 0; i < n; ++i)
        {
            Vq = solveVp (biasVg, (float) cathodeVk);
            cathodeVk = (double) lastIp * cathodeRk;
        }
        outScale = 1.f / std::max (40.f, Bplus - Vq);
        VpPrev = Vq;
    }

    inline float korenIp (float Vg, float Vp, float* dIp_dVp = nullptr) const
    {
        const float S  = std::sqrt (kvb + Vp * Vp);
        const float u  = kp * (1.f / mu + Vg / S);
        // softplus(u) = ln(1+e^u), numerically stable
        const float sp = u > 20.f ? u : std::log1p (std::exp (u));
        const float E1 = Vp / kp * sp;
        const float Ip = E1 > 0.f ? std::pow (E1, ex) / kg : 0.f;

        if (dIp_dVp != nullptr)
        {
            const float sig = u > 20.f ? 1.f : 1.f / (1.f + std::exp (-u));
            // dE1/dVp = sp/kp + (Vp/kp)*sig*(du/dVp), du/dVp = -kp*Vg*Vp/S^3
            const float dE1 = sp / kp - Vg * Vp * Vp * sig / (S * S * S);
            *dIp_dVp = E1 > 0.f ? ex * std::pow (E1, ex - 1.f) / kg * dE1 : 0.f;
        }
        return Ip;
    }

    // Solve f(Vp) = (B+ - Vp)/Rload - Ip(Vg - Vk, Vp - Vk) = 0.
    // Vk is the (slowly varying) cathode voltage; the dVpk/dVp term is 1 so
    // the analytic derivative applies unchanged. Vk == 0 -> legacy solve.
    inline float solveVp (float Vg, float Vk = 0.f)
    {
        float Vp = VpPrev;
        for (int i = 0; i < iterations; ++i)
        {
            float dIp;
            const float Ip = korenIp (Vg - Vk, Vp - Vk, &dIp);
            const float f  = (Bplus - Vp) / Rload - Ip;
            const float df = -1.f / Rload - dIp;   // strictly negative -> safe
            Vp = std::clamp (Vp - f / df, 0.f, Bplus);
            lastIp = Ip;
        }
        return Vp;
    }

    // Returns normalized plate swing: >0 for positive grid drive.
    inline float process (float x)
    {
        const float Vg = biasVg + x * inVolts;
        const float Vp = solveVp (Vg, (float) cathodeVk);
        VpPrev = Vp;
        // Rk||Ck one-pole tracks Ip*Rk -> self-bias sag (off when unconfigured)
        if (cathodeRk > 0.f && cathodePole > 0.0)
            cathodeVk += (1.0 - cathodePole) * ((double) lastIp * cathodeRk - cathodeVk);
        return (Vq - Vp) * outScale;
    }
};

} // namespace eon
