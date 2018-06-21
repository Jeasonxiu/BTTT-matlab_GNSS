function gd=xyz2gd(xyz)
X=xyz(1); Y=xyz(2); Z=xyz(3);

a=6378137.0; f=1/298.257222101;
b=a*(1.-f);

aSq = a^2;
bSq = b^2;
eSq = (aSq - bSq)/aSq;
Lon = atan2(Y,X)*180/pi;

if Lon>180.
    Lon = Lon - 360.;
elseif Lon < -180.
    Lon = Lon+360.;
end

p=sqrt(X^2 + Y^2);
q=0;
Phi0 = atan2(Z*inv(1-eSq),p);

while(q ~= 1)
    N0 = aSq/sqrt(aSq*(cos(Phi0))^2 + bSq*(sin(Phi0))^2);
    h = p/cos(Phi0) - N0;
    Phi = atan2(Z*inv(1-eSq*(N0/(N0+h))),p);
    if abs(Phi - Phi0) <= 1e-13
        break;
    else
        Phi0 = Phi;
    end
end
Lat = Phi * 180/pi;
gd = [Lat Lon h];