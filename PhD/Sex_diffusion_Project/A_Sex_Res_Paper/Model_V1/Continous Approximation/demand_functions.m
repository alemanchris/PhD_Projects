clear all
[sigma,endow,neg,a,na,beta,pp,eh_g,el_g,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,ngk,alpha,nit,tol,epsilon,tol1,gamma,phy]=parameters(1);
index_neg = ag<0;
index_pos = ag>0;
neg_ag    = ag.*index_neg;
pos_ag    = ag.*index_pos;
stream_p = linspace(0.1,0.8,6);
%stream_r = linspace(0.015,0.022,6);
stream_r = linspace(0.015,0.0185,6);
sex_d = zeros(1,size(stream_p,2));
sex_s = zeros(1,size(stream_p,2));
asset_d = zeros(1,size(stream_r,2));
asset_s = zeros(1,size(stream_r,2));

for i = 1:size(stream_p,2)
    %{
    [e_m_e,~]   = excess([0.017,stream_p(i)],1,1);
    [e_m_ne,~]  = excess([0.017,stream_p(i)],2,1);
    [e_f_e,~]   = excess([0.017,stream_p(i)],1,2);
    [e_f_ne,~]  = excess([0.017,stream_p(i)],2,2);
    sex_d(i) = sum(sum([e_m_e(2),e_m_ne(2)]));
    sex_s(i) = sum(sum([e_f_e(2),e_f_ne(2)]));
    %}
    % average asset demand and average asset supply
    [~,gk_m_e]   = excess([stream_r(i),0.4],1,1);
    [~,gk_m_ne]  = excess([stream_r(i),0,4],2,1);
    [~,gk_f_e]   = excess([stream_r(i),0,4],1,2);
    [~,gk_f_ne]  = excess([stream_r(i),0,4],2,2);
    asset_d(i) = abs(neg_ag*sum(gk_m_e,2)+neg_ag*sum(gk_m_ne,2)+neg_ag*sum(gk_f_e,2)+neg_ag*sum(gk_f_ne,2)); % Do not reweight
    asset_s(i) = pos_ag*sum(gk_m_e,2)+pos_ag*sum(gk_m_ne,2)+pos_ag*sum(gk_f_e,2)+pos_ag*sum(gk_f_ne,2);
end
%%
filename = 'curves.mat';
save(filename)
%%
clear all
load('curves.mat')
%%
figure(14)
plot(sex_d,stream_p,'-r',sex_s,stream_p,'-b');
title('Sex market');
xlabel('Extramarital risky sex');
ylabel('Price');
%%
figure(15)
plot(abs(asset_d),stream_r,'-r',asset_s,stream_r,'-b');
title('Asset market');
xlabel('Assets');
ylabel('Interest rate (r)');