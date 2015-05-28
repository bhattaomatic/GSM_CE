function y = gaussian_pulse_shaping_filter(BT,sps,T)
    
    %T = 0.05;
    %sps = 18;
    Ts = T/sps;
    t = -2*T:Ts:2*T;
    B = BT/T;
    alpha = (2*pi*B)/sqrt(log(2));
    gauss = ((2*T)^-1)*(qfunc(alpha*(2*t-(T/2))) - qfunc(alpha*(2*t+(T/2))));
    K = (pi/2)/sum(gauss);
    y = K * gauss;
    %plot(t/T,y);
end
