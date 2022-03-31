function figPTF()
  % this function uses the early pedo-transfer functions of Rawls, W. J.,
  % & Brakensiek, D. L. (1989). Estimation of Soil Water Retention and 
  % Hydraulic Properties. In H. J. Morel-Seytoux (Ed.), Unsaturated Flow in
  % Hydrologic Modeling (pp. 275–300). Springer Netherlands. 
  % https://doi.org/10.1007/978-94-009-2352-2_10
  % and derives even more simple relationships that can be used to infer
  % unsaturated zone parameters (alpha/N for van Genuchten; alternatively
  % air entry pressure and pore size index for Brooks-Corey) from saturated
  % hydraulic conductivity alone. This is of course very simplistic, but
  % might help in parameter/dimension reduction for model calibration.
  
  % we start with creating a "random" sample within ternary coordinates, to
  % evaluate Rawls/Brakensiek across all the parameter space.
  
  % create random numbers and ensure that (X+Y+Z) == 1 (<- ternary)
  nP        = 1500;
  hs        = haltonset(3,'Skip',10*nP,'Leap',nthprime(3+3)-1);
  hs        = scramble(hs,'RR2');
  r         = hs(1:nP,:);
  r         = -log(1-r);
  r         = r./sum(r,2);
  
  % convert ternary coordinates to cartesian ones for plotting
  [rX, rY]  = tern2cart(r(:,1),r(:,2));
  
  % make ternary plot
  pltWdth         = 14.4;
  pltHght         = 7;

  fontname        = 'Cochineal';
  fontsize        = 10.95;
  fig             = figure(1);
  fig             = clf(fig);
  fig.Units       = "centimeters";
  fig.PaperUnits  = "centimeters";
  fig.PaperSize   = [pltWdth,pltHght];
  fig.Position    = [0,0,pltWdth,pltHght];

  fig             = clf(fig);
  set(fig,'DefaultAxesFontSize',fontsize,'DefaultAxesFontName', fontname);

  t               = tiledlayout(1,2);
  t.TileSpacing   = "compact";
  t.Padding       = "compact";

  nexttile([1 1])
  nCells          = 5;
  dL              = 0.08;
  [x, y]          = gridcoords( nCells, dL );
  
  % draw the ternary coordinate system
  line("XData", x,"YData", y,"Color", [0.75 0.75 0.75 1], "LineWidth", 1.5 );
  line( "XData", [-0.5, 0.5, 0, -0.5], ...
        "YData", [0, 0, sqrt( 3 )/2, 0], ...
        "ZData", [0, 0, 0, 0], ...
        "Color", "k", ...
        "LineWidth", 2);
  text(0, -2*dL, "← sand", "HorizontalAlignment", "center" ,...
        "FontSize",fontsize, "Clipping","on",'FontName',fontname);
  text( -0.25 - 2*dL*sqrt(0.75), sqrt( 3 )/4 + dL, "clay →", ...
        "HorizontalAlignment", "center", ...
        "Rotation", 60,"FontSize",fontsize , "Clipping","on",'FontName',fontname);
  text(0.25+2*dL*sqrt(0.75), sqrt( 3 )/4 + dL, "silt →", ...
        "HorizontalAlignment", "center", ...
        "Rotation", -60 ,"FontSize",fontsize, "Clipping","on",'FontName',fontname);
  axis off
  tP = 1:1:nCells+1;
  
  text(1/2 - (tP-1) ./ nCells, ...
    - dL * ones( 1, length( tP ) ), ...
    compose('%.0f', round( 100*(tP-1) ./ nCells, 2 )' ), ...
    "HorizontalAlignment", "left" , "Clipping","on",'fontsize',fontsize, ...
    'FontName',fontname);
  text(-1/2 + 1/2 * (tP-1)./nCells - dL*sqrt(0.75), ...0
    sqrt( 3 ) / 2 * (tP-1) / nCells + dL*0.5, ...
    compose('%.0f', round( 100*(tP-1) ./ nCells, 2 )' ), ...
    "HorizontalAlignment", "right", ...
    "Rotation", 60 , "Clipping","on",'fontsize',fontsize, ...
    'FontName',fontname);
  text(1/2 * (tP-1) ./ nCells + sqrt(0.75)*dL, ...
    sqrt( 3 ) / 2 * (1 - (tP-1) ./ nCells) + dL*0.5, ...
    compose('%.0f', round( 100*(tP-1) ./ nCells, 2 )' ), ...
    "HorizontalAlignment", "right", ...
    "Rotation", -60 , "Clipping","on",'fontsize',fontsize, ...
    'FontName',fontname);
  
  % plot the random sample; it should have equal density everywhere
  hold on
  plot(rX,rY,'.',"Color",[0.6 0.6 0.6])

  % add soil classification info
  pol = {[ 100, 0; 90, 10; 85, 0], "sand";
          [85,0; 90,10; 85,15; 70,0],"loamy sand";
          [70,0; 85,15; 80,20; 52.5,20; 52.5,7.5; 42.5,7.5; 50,0],"sandy loam";
          [50,0; 22.5,27.5; 0,27.5; 0,12.5; 7.5,12.5; 20,0],"silt loam";
          [0,0; 20,0; 7.5,12.5; 0,12.5],"silt"
          [42.5,7.5; 52.5,7.5; 52.5,20; 45,27.5; 22.5,27.5],"loam";
          [52.5,20; 80,20; 65,35; 45,35; 45,27.5 ],"sandy clay loam";
          [20,27.5; 45,27.5; 45,40; 20,40],"clay loam";
          [0,27.5; 20,27.5; 20,40; 0,40],"silty clay loam";
          [0,40; 20,40; 0,60],"silty clay";
          [45,35; 65,35; 45,55],"sandy clay";
          [20,40; 45,40; 45,55; 0,100; 0,60],"clay"};
  
  for j = 1:size(pol,1)
    [xPol, yPol]  = tern2cart( pol{j}(:,1)/100, pol{j}(:,2)/100);
    poly          = polyshape(xPol,yPol);
    plot(poly,'FaceAlpha',0.0,'FaceColor',[0.8 0.8 0.8],'LineWidth',1.5)
  end

  % as this does not use Matlab's default axis, we need some manual changes
  ylim('manual')
  daspect([1 1 1])
  set(gca,'Clipping','off');
  zoom(1.1)

  % define sand and clay contents for rawlsBrakensiek
  S       = r(:,1)*100;
  C       = r(:,2)*100;
  
  % determine porosity from bulk density estimator; apply some randomness
  BD      = 1.51+0.0025*S-0.0013*S*1.5-0.0006*C*1.5;
  P       = (2.65-BD)/2.65;
  hs      = haltonset(1,'Skip',10*nP,'Leap',nthprime(1+3)-1);
  hs      = scramble(hs,'RR2');
  r       = hs(1:nP,:);
  P       = P+0.15*(r-0.5);
  
  % rawlsBrakensiek estimators
  Ks      = +19.52348*P-8.96847-0.028212*C+0.00018107*S.*S...
            -0.0094125*C.*C-8.395215*P.*P+0.077718.*S.*P-0.00298*S.*S.*P.*P...
            -0.019492*C.*C.*P.*P+0.0000173*S.*S.*C+0.02733*C.*C.*P...
            +0.001434*S.*S.*P-0.0000035*C.*C.*S;
  Ks      = exp(Ks)/100/3600;
  
  lam     = -0.7842831+0.0177544*S-1.062498*P-0.00005304*S.*S...
            -0.00273493*C.*C+1.11134946*P.*P-0.03088295*S.*P...
            +0.00026587*S.*S.*P.*P-0.00610522*C.*C.*P.*P...
            -0.00000235*S.*S.*C+0.00798746*C.*C.*P-0.00674491*P.*P.*C;
  lam     = exp(lam);
  
  hB      = +5.3396738+0.1845038*C-2.48394546*P-0.00213853*C.*C...
            -0.04356349*S.*P-0.61745089*C.*P+0.00143598*S.*S.*P.*P...
            -0.00855375*C.*C.*P.*P-0.00001282*S.*S.*C+0.00895359*C.*C.*P...
            -0.00072472*S.*S.*P+0.0000054*C.*C.*S+0.50028060*P.*P.*C;
  hB      = exp(hB)/100;
  
  % determine a relationship between alpha (= 1/hB) and Ks (log&log)
  yM      = log10(1./hB);
  xM      = log10(Ks);
  
  % sort and smooth
  temp    = sortrows([xM,yM]);
  xS      = temp(:,1);
  yS      = temp(:,2);
  yS      = smoothdata(yS,'lowess','SamplePoints',xS);
  
  % fit & print a sigmoid model
  mdl     = @(b,x) b(1)+(b(2)-b(1))*1./( 1+exp( -b(3).*(x-b(4)) ) );
  mdl     = fitnlm(xM,yM,mdl,[min(yS) max(yS) 1 -7]);
  disp(mdl)
  xE      = linspace(-12, -2,100);
  logAlph = mdl.predict(xE');
%   logAlph =  -0.97 +(-0.97-5.00)*1./( 1+exp( -0.34.*(xE'+2.74) ) );
  
  % plot it
  nexttile([1 1])
  scatter(10.^xM,10.^yM,5,[0.6 0.6 0.6],'filled','DisplayName','Rawls & Brakensiek (1989)')
  hold on
%   loglog(10.^xS,10.^yS,'b.','LineWidth',2)
  plot(10.^xE,10.^logAlph,"Color",[0.8 0.1 0.3],'LineWidth',2,'DisplayName','inferred relationship');
  xlim([10.^min(xE) 10.^max(xE)])
  xlabel("{\it K}_{sat}")
  ylabel("{\it α} in 1/m^{-1}")
  set(gca,"FontSize",fontsize)
  set(gca,'YScale','log')
  set(gca,'XScale','log')
  grid on
  box on
  set(gca,"MinorGridLineStyle","none")
  
  leg = legend('Location','best','fontsize',fontsize);
  leg.ItemTokenSize = [12 18];
  leg.EdgeColor     = 'none';
  leg.Color         = 'none';

%   % determine a relationship between N (Lambda+1) and Ks (log&log)
%   xM      = log10(Ks);
%   yM      = log10(lam);
%   
%   % sort and smooth
%   temp    = sortrows([xM,yM]);
%   xS      = temp(:,1);
%   yS      = temp(:,2);
%   yS      = smoothdata(yS,'lowess','SamplePoints',xS);
%   
%   % fit & print a sigmoid model
%   mdl     = @(b,x) b(1)+(b(2)-b(1))*1./( 1+exp( -b(3).*(x-b(4)) ) );
%   mdl     = fitnlm(xM,yM,mdl,[min(yS) max(yS) 1 -7]);
%   disp(mdl)
%   logLam  = mdl.predict(xE');
  
%   % plot it
%   nexttile([2 4])
%   semilogx(10.^xM,10.^yM+1,'.',"Color",[0.5 0.5 0.5 0.75])
%   hold on
% %   semilogx(10.^xS,10.^yS+1,'b.','LineWidth',2)
%   semilogx(10.^xE,10.^logLam+1,"Color",[0.8 0.1 0.3],'LineWidth',2);
%   xlim([10.^min(xE) 10.^max(xE)])
%   xlabel("{\it K}_{sat}")
%   ylabel("{\it N}")
%   set(gca,"FontSize",14)
%   grid on
%   set(gca,"MinorGridLineStyle","none")

  % make sure the full plot is saved
  axT = axes;
  axT.Position = [0,0,1,1];
  axT.Color = 'None';%
  axT.XTick = [0 1];
  axT.XColor = 'white';
  axT.YTick = [0 1];
  axT.YColor = 'white';

  exportgraphics(fig,"./figure.pdf",'ContentType','vector')
end

function [x, y] = gridcoords( res,dL )
  %GRIDCOORDS Compute the grid coordinates.
  
  % Create an equally-spaced vector.
  t = linspace(0,1,res+1);

  % First, assemble the x-coordinates.
  x1 = [ 1/2 - t + 0.5*dL*sqrt(0.25);
        -1/2 * t;
        NaN( 1, res+1 )];
  x2 = [1/2 - t;
    1/2 - 1/2 * t + 0.5*dL*sqrt(0.25);
    NaN( 1, res+1 )];
  x3 = [-1/2 + 1/2 * t - 0.5*dL;
    1/2 - 1/2 * t;
    NaN( 1, res+1 )];
  x = [x1; x2; x3];
  x = x(:);
  % Next, assemble the y-coordinates.
  y1 = [zeros( 1, res+1 )- 0.5*dL*sqrt(0.75);
    sqrt( 3 ) / 2 * (1 - t);
    NaN( 1, res+1 )];
  y2 = [zeros( 1, res+1 );
    sqrt( 3 ) / 2 * t + 0.5*dL*sqrt(0.75);
    NaN( 1, res+1 )];
  y3 = [sqrt( 3 ) / 2 * t ;
    sqrt( 3 ) / 2 * t;
    NaN( 1, res+1 )];
  y = [y1; y2; y3];
  y = y(:);
end

function [x, y] = tern2cart( A, B)
  %TERN2CART Convert the ternary coordinates (A, B, C) to
  %Cartesian coordinates (x, y).
  
  y = B * sind( 60 );
  x = 1 - A - y * cotd( 60 ) - 1/2;  
end

