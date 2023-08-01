#ifndef WDF_H_INCLUDED
#define WDF_H_INCLUDED

#include "omega.h"
#include <string>
#include <cmath>

namespace WaveDigitalFilter
{
/** Wave digital filter base class */
template <typename T>
class WDF
{
public:
    WDF (std::string type) : type (type) {}
    virtual ~WDF() {}

    virtual inline void calcImpedance() {}
    virtual inline void propagateImpedance() {}

    virtual inline void incident (T x) noexcept {}
    virtual inline T reflected() noexcept { return 0.0; }


    inline T voltage() const noexcept
    {
        return (a + b) / 2.0;
    }

    inline T current() const noexcept
    {
        return (a - b) / (2.0 * R);
    }

public:
    T a = 0.0; // incident wave
    T b = 0.0; // reflected wave
    T R = 1.0e-9;
    T G = 1.0 / R;

private:
    const std::string type;
};

/** WDF node base class */
template <typename T>
class WDFNode : public WDF<T>
{
public:
    WDFNode (std::string type) : WDF<T> (type) {}
    virtual ~WDFNode() {}

    void connectToNode (WDF<T>* node)
    {
        next = node;
    }

    inline void propagateImpedance() override
    {
        calcImpedance();

        if (next != nullptr)
            next->propagateImpedance();
    }

protected:
    WDF<T>* next = nullptr;
};

/** WDF Resistor Node */
template <typename T>
class Resistor : public WDFNode<T>
{
public:
    Resistor (T value) :
        WDFNode<T> ("Resistor"),
        R_value (value)
    {
        calcImpedance();
    }
    virtual ~Resistor() {}

    void setResistanceValue (T newR)
    {
        if (newR == R_value)
            return;

        R_value = newR;
        propagateImpedance();
    }

    inline void calcImpedance() override
    {
        R = R_value;
        G = 1.0 / R;
    }

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = 0.0;
        return b;
    }

private:
    T R_value = 1.0e-9;
};

/** WDF Capacitor Node */
template <typename T>
class Capacitor : public WDFNode<T>
{
public:
    Capacitor (T value, T fs, T alpha = 1.0) :
        WDFNode<T> ("Capacitor"),
        C_value (value),
        fs (fs),
        alpha (alpha),
        b_coef ((1.0 - alpha) / 2.0),
        a_coef ((1.0 + alpha) / 2.0)
    {
        calcImpedance();
    }
    virtual ~Capacitor() {}

    void setCapacitanceValue (T newC)
    {
        if (newC == C_value)
            return;

        C_value = newC;
        propagateImpedance();
    }

    inline void calcImpedance() override
    {
        R = 1.0 / ((1.0 + alpha) * C_value * fs);
        G = 1.0 / R;
    }

    inline void incident (T x) noexcept override
    {
        a = x;
        z = a;
    }

    inline T reflected() noexcept override
    {
        b = b_coef * b + a_coef * z;
        return b;
    }

private:
    T C_value = 1.0e-6;
    T z = 0.0;

    const T fs;
    const T alpha;

    const T b_coef;
    const T a_coef;
};


/** WDF Inductor Node */
template <typename T>
class Inductor : public WDFNode<T>
{
public:
    Inductor (T value, T fs, T alpha = 1.0) :
        WDFNode<T> ("Inductor"),
        L_value (value),
        fs (fs),
        alpha (alpha),
        b_coef ((1.0 - alpha) / 2.0),
        a_coef ((1.0 + alpha) / 2.0)
    {
        calcImpedance();
    }
    virtual ~Inductor() {}

    void setInductanceValue (T newL)
    {
        if (newL == L_value)
            return;

        L_value = newL;
        propagateImpedance();
    }

    inline void calcImpedance() override
    {
        R = (1.0 + alpha) * L_value * fs;
        G = 1.0 / R;
    }

    inline void incident (T x) noexcept override
    {
        a = x;
        z = a;
    }

    inline T reflected() noexcept override
    {
        b = b_coef * b - a_coef * z;
        return b;
    }

private:
    T L_value = 1.0e-6;
    T z = 0.0;

    const T fs;
    const T alpha;

    const T b_coef;
    const T a_coef;
};

/** WDF Switch */
template <typename T>
class Switch : public WDFNode<T>
{
public:
    Switch():
        WDFNode ("Switch")
    {}
    virtual ~Switch() {}

    inline void calcImpedance() override {}

    void setClosed (bool shouldClose) { closed = shouldClose; }

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = closed ? -a : a;
        return b;
    }

private:
    bool closed = true;
};

/** WDF Open */
template <typename T>
class Open : public WDFNode<T>
{
public:
    Open():
        WDFNode ("Open")
    {}
    virtual ~Open()
    {
        R = 1.0e15;
        G = 1.0 / R;
    }

    inline void calcImpedance() override {}

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = a;
        return b;
    }
};

/** WDF Short */
template <typename T>
class Short : public WDFNode<T>
{
public:
    Short():
        WDFNode<T> ("Short")
    {}
    virtual ~Short()
    {
        R = 1.0e-15;
        G = 1.0 / R;
    }

    inline void calcImpedance() override {}

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = -a;
        return b;
    }
};

/** WDF Voltage Polarity Inverter */
template <typename T>
class PolarityInverter : public WDFNode<T>
{
public:
    PolarityInverter (WDFNode* port1) :
        WDFNode<T> ("Polarity Inverter"),
        port1 (port1)
    {
        port1->connectToNode (this);
        calcImpedance();
    }
    virtual ~PolarityInverter() {}

    inline void calcImpedance() override
    {
        R = port1->R;
        G = 1.0 / R;
    }

    inline void incident (T x) noexcept override
    {
        a = x;
        port1->incident (-x);
    }

    inline T reflected() noexcept override
    {
        b = -port1->reflected();
        return b;
    }

private:
    WDFNode* port1;
};

/** WDF y-parameter 2-port (short circuit admittance) */
template <typename T>
class YParameter : public WDFNode<T>
{
public:
    YParameter (WDFNode* port1, T y11, T y12, T y21, T y22) :
        WDFNode<T> ("YParameter"),
        port1 (port1)
    {
        y[0][0] = y11; y[0][1] = y12;
        y[1][0] = y21; y[1][1] = y22;

        port1->connectToNode (this);
        calcImpedance();
    }

    virtual ~YParameter() {}

    inline void calcImpedance() override
    {
        denominator = y[1][1] + port1->R * y[0][0] * y[1][1] - port1->R * y[0][1] * y[1][0];
        R = (port1->R * y[0][0] + 1.0) / denominator;
        G = 1.0 / R;

        T rSq = port1->R * port1->R;
        T num1A = -y[1][1] * rSq * y[0][0] * y[0][0];
        T num2A = y[0][1] * y[1][0] * rSq * y[0][0];

        A = (num1A + num2A + y[1][1]) / (denominator * (port1->R * y[0][0] + 1.0));
        B = -port1->R * y[0][1] / (port1->R * y[0][0] + 1.0);
        C = -y[1][0] / denominator;
    }

    inline void incident (T x) noexcept override
    {
        a = x;
        port1->incident(A * port1->b + B * x);
    }

    inline T reflected() noexcept override
    {
        b = C * port1->reflected();
        return b;
    }

private:
    WDFNode* port1;
    T y[2][2] = {{ 0.0, 0.0 }, { 0.0, 0.0 }};

    T denominator = 1.0;
    T A = 1.0f;
    T B = 1.0f;
    T C = 1.0f;
};

/** WDF 3-port adapter base class */
template <typename T>
class WDFAdaptor : public WDFNode<T>
{
public:
    WDFAdaptor (WDFNode<T>* port1, WDFNode<T>* port2, std::string type) :
        WDFNode<T> (type),
        port1 (port1),
        port2 (port2)
    {
        port1->connectToNode (this);
        port2->connectToNode (this);
    }
    virtual ~WDFAdaptor() {}

protected:
    WDFNode* port1;
    WDFNode* port2;
};

/** WDF 3-port parallel adaptor */
template <typename T>
class WDFParallel : public WDFAdaptor<T>
{
public:
    WDFParallel (WDFNode<T>* port1, WDFNode<T>* port2) :
        WDFAdaptor<T> (port1, port2, "Parallel")
    {
        calcImpedance();
    }
    virtual ~WDFParallel() {}

    inline void calcImpedance()
    {
        G = port1->G + port2->G;
        R = 1.0 / G;
        port1Reflect = port1->G/G;
        port2Reflect = port2->G/G;
    }

    inline T reflected() noexcept override
    {
        b = port1Reflect * port1->reflected() + port2Reflect * port2->reflected();
        return b;
    }

    inline void incident (T x) noexcept override
    {
        port1->incident (x + (port2->b - port1->b) * port2Reflect);
        port2->incident (x + (port2->b - port1->b) * -port1Reflect);
        a = x;
    }

private:
    T port1Reflect = 1.0;
    T port2Reflect = 1.0;
};

/** WDF 3-port series adaptor */
template <typename T>
class WDFSeries : public WDFAdaptor<T>
{
public:
    WDFSeries (WDFNode<T>* port1, WDFNode<T>* port2) :
        WDFAdaptor<T> (port1, port2, "Series")
    {
        calcImpedance();
    }
    virtual ~WDFSeries() {}

    inline void calcImpedance()
    {
        R = port1->R + port2->R;
        G = 1.0 / R;
        port1Reflect = port1->R/R;
        port2Reflect = port2->R/R;
    }

    inline T reflected() noexcept override
    {
        b = -(port1->reflected() + port2->reflected());
        return b;
    }

    inline void incident (T x) noexcept override
    {
        port1->incident (port1->b - port1Reflect * (x + port1->b + port2->b));
        port2->incident (port2->b - port2Reflect * (x + port1->b + port2->b));

        a = x;
    }

private:
    T port1Reflect = 1.0;
    T port2Reflect = 1.0;
};

/** WDF Voltage source with resistance */
template <typename T>
class ResistiveVoltageSource : public WDFNode<T>
{
public:
    ResistiveVoltageSource (T value = 1.0e-9) :
        WDFNode<T> ("Resistive Voltage"),
        R_value (value)
    {
        calcImpedance();
    }
    virtual ~ResistiveVoltageSource() {}

    void setResistanceValue (T newR)
    {
        if (newR == R_value)
            return;

        R_value = newR;
        propagateImpedance();
    }

    inline void calcImpedance()
    {
        R = R_value;
        G = 1.0 / R;
    }

    void setVoltage (T newV) { Vs = newV; }

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = Vs;
        return b;
    }

private:
    T Vs;
    T R_value = 1.0e-9;
};

/** WDF Voltage source with 1 pOhm resistance */
template <typename T>
class IdealVoltageSource : public WDFNode<T>
{
public:
    IdealVoltageSource() : WDFNode<T> ("IdealVoltage")
    {
        calcImpedance();
    }
    virtual ~IdealVoltageSource() {}

    inline void calcImpedance() {}

    void setVoltage (T newV) { Vs = newV; }

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = -a + 2.0 * Vs;
        return b;
    }

private:
    T Vs;
};

/** WDF Current source with resistance */
template <typename T>
class ResistiveCurrentSource : public WDFNode<T>
{
public:
    ResistiveCurrentSource (T value=1.0e9) :
        WDFNode<T> ("Resistive Current"),
        R_value (value)
    {
        calcImpedance();
    }
    virtual ~ResistiveCurrentSource() {}

    void setResistanceValue (T newR)
    {
        if (newR == R_value)
            return;

        R_value = newR;
        propagateImpedance();
    }

    inline void calcImpedance()
    {
        R = R_value;
        G = 1.0 / R;
    }

    void setCurrent (T newI) { Is = newI; }

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = 2 * R * Is;
        return b;
    }

private:
    T Is;
    T R_value = 1.0e9;
};

/** WDF Current source with 1 GOhm resistance */
template <typename T>
class IdealCurrentSource : public WDFNode<T>
{
public:
    IdealCurrentSource() : WDFNode<T> ("Ideal Current")
    {
        calcImpedance();
    }
    virtual ~IdealCurrentSource() {}

    inline void calcImpedance()
    {
        R = 1.0e9;
        G = 1.0 / R;
    }

    void setCurrent (T newI) { Is = newI; }

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        b = 2 * next->R * Is + a;
        return b;
    }

private:
    T Is;
};

template <typename T> inline int signum (T val)
{
    return (T (0) < val) - (val < T (0));
}

/**
 * Diode pair nonlinearity evaluated in the wave domain
 * See Werner et al., "An Improved and Generalized Diode Clipper Model for Wave Digital Filters"
 * https://www.researchgate.net/publication/299514713_An_Improved_and_Generalized_Diode_Clipper_Model_for_Wave_Digital_Filters
 * */
template <typename T>
class DiodePair : public WDFNode<T>
{
public:
    DiodePair (T Is, T Vt) :
        WDFNode<T> ("DiodePair"),
        Is (Is),
        Vt (Vt)
    {}

    virtual ~DiodePair() {}

    inline void calcImpedance() {}

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        // See eqn (18) from reference paper
        T lambda = (T) signum (a);
        b = a + 2 * lambda * (next->R * Is - Vt * omega4 (float (log (next->R * Is / Vt) + (lambda * a + next->R * Is) / Vt)));
        return b;
    }

private:
    const T Is; // reverse saturation current
    const T Vt; // thermal voltage
};

template <typename T>
class Diode : public WDFNode<T>
{
public:
    Diode (T Is, T Vt) :
        WDFNode<T> ("Diode"),
        Is (Is),
        Vt (Vt)
    {}

    virtual ~Diode() {}

    inline void calcImpedance() {}

    inline void incident (T x) noexcept override
    {
        a = x;
    }

    inline T reflected() noexcept override
    {
        // See eqn (10) from reference paper
        b = a + 2 * next->R * Is - 2 * Vt * omega4 (float (log (next->R * Is / Vt) + (a + next->R * Is) / Vt));
        return b;
    }

private:
    const T Is; // reverse saturation current
    const T Vt; // thermal voltage
};

}

#endif // WDF_H_INCLUDED
