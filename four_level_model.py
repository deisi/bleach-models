import numpy as np
from scipy.special import erf, erfc


def trace(t, Amp, t1, t2, c, mu, sigma):
    """Bleach trace, derived from an analytical Solution of the four level system.

    Function is derived by analytically solving the 4 level system and
    subsequent convolution with a gaussian excitation function of the model.
    Initial state is N0=1, N1=0, N2=0 and N3=0.

    Note:
        This implementation has a problem when t1=t2. Due to numerical constrains
        this must be avoided, as it divides by 0.

    Args:
        t: Array of Time values. Defines units of ``t1`` and ``t2``.
        Amp: Amplitude of the excitation pulse. Determines the fraction
            of oscillators excited by the excitation pulse into N1 state.
        t1: Lifetime of the first excited vibrational state in units of
              ``t``
        t2: Lifetime of the second excited vibrational state in units of
            ``t``.
        c: Scaling factor of final (heated) state. Used to account for
            spectral differences induced by residual heat.
        mu: Tempoaral position of pump pulse in units of ``t``.
        sigma: Temporal width of pump pulse in units of ``t``.

    Returns:
        Modeled result as deduced from the 4 level system for the given array of
        ``t`` time values.

    """
    pi=np.pi;
    #a0 = erf((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/sigma))
    def mysqrt(x): return np.sqrt(x)
    aux0=sigma*((t1**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma))))));
    aux1=sigma*((t1**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma))))));
    aux2=sigma*((t1**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma))))));
    aux3=(((t1-t2)**2))*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma)))));
    aux4=sigma*(t1*(t2*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma)))))));
    aux5=sigma*(t1*(t2*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma)))))));
    aux6=sigma*(t1*(t2*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma)))))));
    aux7=sigma*((t2**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma))))));
    aux8=sigma*((t2**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma))))));
    aux9=sigma*((t2**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)/\
    sigma))))));
    aux10=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux11=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*((t1**2)*(-1.+(erf(aux10))))));
    aux12=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux13=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux12)))))));
    aux14=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux15=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux14)))))));
    aux16=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux17=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux16))))));
    aux18=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux19=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux18))))));
    aux20=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t1))-(((2.**-0.5)\
    *t)/sigma);
    aux21=(np.exp(((2.*((sigma**2)*(t1**-2.)))+(((2.*mu)/t1)+((-2.*t)/t1))\
       )))*((mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux20)))))));
    aux22=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t1))-(((2.**-0.5)\
    *t)/sigma);
    aux23=(np.exp(((2.*((sigma**2)*(t1**-2.)))+(((2.*mu)/t1)+((-2.*t)/t1))\
       )))*((mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux22)))))));
    aux24=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t1))-(((2.**-0.5)\
    *t)/sigma);
    aux25=(np.exp(((2.*((sigma**2)*(t1**-2.)))+(((2.*mu)/t1)+((-2.*t)/t1))\
       )))*((mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux24))))));
    aux26=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t2))-(((2.**-0.5)*\
    t)/sigma);
    aux27=(np.exp((((0.5*((sigma**2)*(t2**-2.)))+(mu/t2))-(t/t2))))*((\
    mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux26)))))));
    aux28=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t2))-(((2.**-0.5)*\
    t)/sigma);
    aux29=(np.exp((((0.5*((sigma**2)*(t2**-2.)))+(mu/t2))-(t/t2))))*((\
    mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux28))))));
    aux30=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t2))-(((2.**-0.5)*\
    t)/sigma);
    aux31=(np.exp((((0.5*((sigma**2)*(t2**-2.)))+(mu/t2))-(t/t2))))*((\
    mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux30))))));
    aux32=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t2))-(((2.**-0.5)*\
    t)/sigma);
    aux33=(np.exp((((0.5*((sigma**2)*(t2**-2.)))+(mu/t2))-(t/t2))))*((\
    mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux32))))));
    aux34=(0.5*((sigma**2)*(t1**-2.)))+((mu/t1)+((0.5*((sigma**2)*(t2**-2.\
       )))+((mu/t2)+(((sigma**2)/t2)/t1))));
    aux35=(((2.**-0.5)*mu)/sigma)+((((2.**-0.5)*sigma)/t1)+(((2.**-0.5)*\
    sigma)/t2));
    aux36=(mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf((aux35-(((2.**-0.5)*t)/sigma))))))));
    aux37=(0.5*((sigma**2)*(t1**-2.)))+((mu/t1)+((0.5*((sigma**2)*(t2**-2.\
       )))+((mu/t2)+(((sigma**2)/t2)/t1))));
    aux38=(((2.**-0.5)*mu)/sigma)+((((2.**-0.5)*sigma)/t1)+(((2.**-0.5)*\
    sigma)/t2));
    aux39=(mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf((aux38-(((2.**-0.5)*t)/sigma)))))));
    aux40=(0.5*((sigma**2)*(t1**-2.)))+((mu/t1)+((0.5*((sigma**2)*(t2**-2.\
       )))+((mu/t2)+(((sigma**2)/t2)/t1))));
    aux41=(((2.**-0.5)*mu)/sigma)+((((2.**-0.5)*sigma)/t1)+(((2.**-0.5)*\
    sigma)/t2));
    aux42=(mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf((aux41-(((2.**-0.5)*t)/sigma)))))));
    aux43=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t2))-(((2.**-0.5)\
    *t)/sigma);
    aux44=(np.exp(((2.*((sigma**2)*(t2**-2.)))+(((2.*mu)/t2)+((-2.*t)/t2))\
       )))*((mysqrt((2.*pi)))*(sigma*((t2**2)*(-1.+(erf(aux43))))));
    aux45=t1*(t2*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t1))-(t*t1)))/t1)/\
    sigma))));
    aux46=(np.exp((0.5*((t1**-2.)*((sigma**2)+((2.*(mu*t1))+(-2.*(t*t1))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux45));
    aux47=t1*(t2*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t1))-(t*t1)))/t1)/\
    sigma))));
    aux48=(np.exp((0.5*((t1**-2.)*((sigma**2)+((2.*(mu*t1))+(-2.*(t*t1))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux47));
    aux49=(t2**2)*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t1))-(t*t1)))/t1)/\
    sigma)));
    aux50=(np.exp((0.5*((t1**-2.)*((sigma**2)+((2.*(mu*t1))+(-2.*(t*t1))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux49));
    aux51=t1*(t2*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t2))-(t*t2)))/t2)/\
    sigma))));
    aux52=(np.exp((0.5*((t2**-2.)*((sigma**2)+((2.*(mu*t2))+(-2.*(t*t2))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux51));
    aux53=(t2**2)*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t2))-(t*t2)))/t2)/\
    sigma)));
    aux54=(np.exp((0.5*((t2**-2.)*((sigma**2)+((2.*(mu*t2))+(-2.*(t*t2))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux53));
    aux55=(3.*(Amp*aux46))+((Amp*(c*aux48))+((-2.*(Amp*aux50))+((Amp*(c*\
    aux52))+(Amp*aux54))));
    aux56=(-2.*((Amp**2)*(c*((np.exp(((aux40-(t/t2))-(t/t1))))*aux42))))+(\
    ((Amp**2)*(c*aux44))+aux55);
    aux57=((Amp**2)*((c**2)*((np.exp(((aux34-(t/t2))-(t/t1))))*aux36)))+((\
    2.*((Amp**2)*((np.exp(((aux37-(t/t2))-(t/t1))))*aux39)))+aux56);
    aux58=((Amp**2)*aux29)+((-2.*((Amp**2)*(c*aux31)))+(((Amp**2)*((c**2)*\
    aux33))+aux57));
    aux59=(2.*((Amp**2)*(c*aux23)))+((-2.*((Amp**2)*aux25))+((2.*((Amp**2)\
    *(c*aux27)))+aux58));
    aux60=(-2.*((Amp**2)*aux17))+((2.*((Amp**2)*(c*aux19)))+((2.*((Amp**2)\
    *aux21))+aux59));
    aux61=((Amp**2)*((c**2)*aux11))+((3.*((Amp**2)*aux13))+((-2.*((Amp**2)\
    *(c*aux15)))+aux60));
    aux62=((Amp**2)*((c**2)*((mysqrt((0.5*pi)))*aux8)))+((Amp*(c*((\
    mysqrt((2.*pi)))*aux9)))+aux61);
    aux63=(2.*((Amp**2)*(c*((mysqrt((2.*pi)))*aux6))))+(((Amp**2)*((\
    mysqrt((0.5*pi)))*aux7))+aux62);
    aux64=(2.*(Amp*((mysqrt((2.*pi)))*aux4)))+((-2.*(Amp*(c*((mysqrt((\
    2.*pi)))*aux5))))+aux63);
    aux65=(Amp*(c*((mysqrt((2.*pi)))*aux2)))+(((mysqrt((0.5*pi)))*(\
    sigma*aux3))+aux64);
    aux66=((Amp**2)*((mysqrt((0.5*pi)))*aux0))+(((Amp**2)*((c**2)*((\
    mysqrt((0.5*pi)))*aux1)))+aux65);
    aux67=(t2**2)*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t2))-(t*t2)))/t2)/\
    sigma)));
    aux68=(np.exp((0.5*((t2**-2.)*((sigma**2)+((2.*(mu*t2))+(-2.*(t*t2))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux67));
    aux69=t1*(t2*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t2))-(t*t2)))/t2)/\
    sigma))));
    aux70=(np.exp((0.5*((t2**-2.)*((sigma**2)+((2.*(mu*t2))+(-2.*(t*t2))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux69));
    aux71=(t1**2)*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t1))-(t*t1)))/t1)/\
    sigma)));
    aux72=(np.exp((0.5*((t1**-2.)*((sigma**2)+((2.*(mu*t1))+(-2.*(t*t1))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux71));
    aux73=(t1**2)*(erfc(((((2.**-0.5)*(((sigma**2)+(mu*t1))-(t*t1)))/t1)/\
    sigma)));
    aux74=(np.exp((0.5*((t1**-2.)*((sigma**2)+((2.*(mu*t1))+(-2.*(t*t1))))\
       ))))*((mysqrt((2.*pi)))*(sigma*aux73));
    aux75=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t2))-(((2.**-0.5)\
    *t)/sigma);
    aux76=(np.exp(((2.*((sigma**2)*(t2**-2.)))+(((2.*mu)/t2)+((-2.*t)/t2))\
       )))*((mysqrt((0.5*pi)))*(sigma*((t2**2)*(-1.+(erf(aux75))))));
    aux77=((((aux66-(Amp*(c*aux68)))-(Amp*aux70))-(Amp*(c*aux72)))-(Amp*\
    aux74))-((Amp**2)*((c**2)*aux76));
    aux78=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t2))-(((2.**-0.5)\
    *t)/sigma);
    aux79=(np.exp(((2.*((sigma**2)*(t2**-2.)))+(((2.*mu)/t2)+((-2.*t)/t2))\
       )))*((mysqrt((0.5*pi)))*(sigma*((t2**2)*(-1.+(erf(aux78))))));
    aux80=(0.5*((sigma**2)*(t1**-2.)))+((mu/t1)+((0.5*((sigma**2)*(t2**-2.)))+((mu/t2)+(((sigma**2)/t2)/t1))));
    aux81=(((2.**-0.5)*mu)/sigma)+((((2.**-0.5)*sigma)/t1)+(((2.**-0.5)*\
    sigma)/t2));
    aux82=(mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf((aux81-(((2.**-0.5)*t)/sigma))))))));
    aux83=(aux77-((Amp**2)*aux79))-((Amp**2)*((np.exp(((aux80-(t/t2))-(t/\
    t1))))*aux82));
    aux84=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t2))-(((2.**-0.5)*\
    t)/sigma);
    aux85=(np.exp((((0.5*((sigma**2)*(t2**-2.)))+(mu/t2))-(t/t2))))*((\
    mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux84)))))));
    aux86=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t2))-(((2.**-0.5)*\
    t)/sigma);
    aux87=(np.exp((((0.5*((sigma**2)*(t2**-2.)))+(mu/t2))-(t/t2))))*((\
    mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux86)))))));
    aux88=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t1))-(((2.**-0.5)\
    *t)/sigma);
    aux89=(np.exp(((2.*((sigma**2)*(t1**-2.)))+(((2.*mu)/t1)+((-2.*t)/t1))\
       )))*((mysqrt((2.*pi)))*(sigma*((t1**2)*(-1.+(erf(aux88))))));
    aux90=((aux83-((Amp**2)*((c**2)*aux85)))-((Amp**2)*aux87))-((Amp**2)*(\
    c*aux89));
    aux91=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t1))-(((2.**-0.5)\
    *t)/sigma);
    aux92=(np.exp(((2.*((sigma**2)*(t1**-2.)))+(((2.*mu)/t1)+((-2.*t)/t1))\
       )))*((mysqrt((0.5*pi)))*(sigma*((t1**2)*(-1.+(erf(aux91))))));
    aux93=((((2.**-0.5)*mu)/sigma)+(((mysqrt(2.))*sigma)/t1))-(((2.**-0.5)\
    *t)/sigma);
    aux94=(np.exp(((2.*((sigma**2)*(t1**-2.)))+(((2.*mu)/t1)+((-2.*t)/t1))\
       )))*((mysqrt((0.5*pi)))*(sigma*((t1**2)*(-1.+(erf(aux93))))));
    aux95=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux96=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*(t1*(t2*(-1.+(erf(aux95)))))));
    aux97=((aux90-((Amp**2)*((c**2)*aux92)))-((Amp**2)*aux94))-((Amp**2)*(\
    (c**2)*aux96));
    aux98=((((2.**-0.5)*mu)/sigma)+(((2.**-0.5)*sigma)/t1))-(((2.**-0.5)*\
    t)/sigma);
    aux99=(np.exp((((0.5*((sigma**2)*(t1**-2.)))+(mu/t1))-(t/t1))))*((\
    mysqrt((2.*pi)))*(sigma*((t1**2)*(-1.+(erf(aux98))))));
    aux100=sigma*((t2**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*\
    t)/sigma))))));
    aux101=sigma*((t2**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*\
    t)/sigma))))));
    aux102=((aux97-((Amp**2)*aux99))-((Amp**2)*(c*((mysqrt((2.*pi)))*\
    aux100))))-(Amp*((mysqrt((2.*pi)))*aux101));
    aux103=sigma*(t1*(t2*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)\
    /sigma)))))));
    aux104=sigma*(t1*(t2*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*t)\
    /sigma)))))));
    aux105=(aux102-((Amp**2)*((c**2)*((mysqrt((2.*pi)))*aux103))))-((\
    Amp**2)*((mysqrt((2.*pi)))*aux104));
    aux106=sigma*((t1**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*\
    t)/sigma))))));
    aux107=sigma*((t1**2)*(1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*\
    t)/sigma))))));
    aux108=(aux105-((Amp**2)*(c*((mysqrt((2.*pi)))*aux106))))-(Amp*((\
    mysqrt((2.*pi)))*aux107));
    aux109=(((t1-t2)**2))*(-1.-(erf(((((2.**-0.5)*mu)/sigma)-(((2.**-0.5)*\
    t)/sigma)))));
    aux110=((2.*pi)**-0.5)*(((t1-t2)**-2.)*(aux108-((mysqrt((0.5*pi)\
       ))*(sigma*aux109))));
    output=aux110/sigma;
    return output
