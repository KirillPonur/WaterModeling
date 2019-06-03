function Kmax=KL_LIMIT(U10,LL)
disp('calling KL_LIMIT function')
NLS= round(1000*LL);
switch NLS
    case  8 
        Kmax=5.95638868922385001+...
            0.0720347041917690456*U10+...
            584.172993830735945/sqrt(U10)-...
            495.241636423855795*exp(-U10);
    case  21
        Kmax=171.844249137808972-...
            29.9715115049259685*U10+...
            11.0462317475020238*U10*log(U10)-...
            0.872695543932866229*U10^2+...
            0.139324811425696058*U10^2*log(U10);
    case  30
        Kmax=111.815140123056191-...
            17.9371868277202707*U10+...
            3.09504183073677383*U10^2-...
            1.2655962097427494*U10^2*log(U10)+...
            0.314107775309787268*U10^2.5;
    case  55
        Kmax=7.49318039534745137-...
            0.1074769345744707*U10+...
            0.494588232199318583*sqrt(U10)+...
            81.5368708904397176*log(U10)/U10+...
            43.8821633213996845/(U10^2);
    case  100
        Kmax=3.50331310640510917-...
            0.13811269308311842*U10+...
            0.00839835244584838482*U10^2-...
            0.000303039793256592249*U10^3+...
            4.47267177226626458e-06*U10^4;
        Kmax= exp(Kmax);
    case   230
        Kmax=2.90160061283935079-...
            0.187283073237593592*sqrt(U10)+...
            17.8301004517745961*log(U10)/U10+...
            10.337585052957506/U10^2;
    otherwise
        error(' SUCH WAVELENGHT IS ABSENT IN PROGRAM ');
end
