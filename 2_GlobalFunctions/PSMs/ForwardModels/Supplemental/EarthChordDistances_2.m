function Chord_Dists_Mat=EarthChordDistances_2(llPoints1,llPoints2)
% Calculates the chordal distance between a set of points on Earth.
%
% INPUTS:
% Each input is N/M by 2, each row a (long, lat) pair, -180<long<180; -90<lat<90.
%
% OUTPUTS:
% Output is a N by M matrix of chordal distances in km (approximating the
% earth as a sphere), where the (i,j) entry is the distance between the ith
% row of llPoints and the jth row of llPoints2.
%
% NOTE: the radius of the earth is taken as 6378.137
%
% function written by Martin P. Tingley, 2012
%

RR=6378.137; %radius of the earth in km
N=length(llPoints1(:,1));
M=length(llPoints2(:,1));

%make a N*M by 4 matrix. Each row one of the N*M possible sets of
%two points:

Pts_paired_Vec=[kron(llPoints1, ones(M,1)), kron(ones(N,1), llPoints2)];
Half_Angles_AsVec=asin(sqrt(sin((Pts_paired_Vec(:,2)-Pts_paired_Vec(:,4))*pi/180/2).^2 + cos(Pts_paired_Vec(:,2)*pi/180).*cos(Pts_paired_Vec(:,4)*pi/180).*sin(abs(Pts_paired_Vec(:,1)-Pts_paired_Vec(:,3))*pi/180/2).^2));
Chords_as_vec=2*RR*sin(Half_Angles_AsVec);

Chord_Dists_Mat=reshape(Chords_as_vec, M, N)';

