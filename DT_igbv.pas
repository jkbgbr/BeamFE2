unit DT_igbv;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  ExtCtrls, Menus, ComCtrls, StdCtrls, Buttons, Tabnotbk, Grids,
  u_mainf, u_2dview, u_igbv_types,u_FormProgress;


//const

procedure StatDocShort(nlc,nbeam,nnode:integer;npoints:intbeamtype;xrel:realbeampointstype;
                       beamtop:beamtopoltype;alfv:realbeamtype;name:shortstring;
                       Nel,Tel,Mel,defx,defy:realbeamlcpointstype);
procedure StatDocFull(nlc,nbeam,nnode:integer;npoints:intbeamtype;xrel:realbeampointstype;
                      beamtop:beamtopoltype;alfv:realbeamtype;name:shortstring;
                      Nel,Tel,Mel,defx,defy:realbeamlcpointstype);
procedure CalcStructure(var beamtop:beamtopoltype;var Av,Iv:realbeamtype;
                        var npoints:intbeamtype;var xrel:realbeampointstype;
                        var Nel,Tel,Mel,defx,defy:realbeamlcpointstype;
                        var reac:reallcnodeproptype;var activecp:beamtopoltype);


implementation

uses
u_formcalculate;


//*************************************************
//            *************************
//**************************************************
//
//            Code-independent Routines
//
//*************************************************
//            *************************
//**************************************************





//*******************
//
//General Subroutines
//
//*******************

function Sign(a:real):integer;

begin
   if a>0 then Sign:=1;
   if a<0 then Sign:=-1;
   if a=0 then Sign:=0;
end; //Sign


function SignLim(a,lim:real):integer;

begin
   if (a>lim) and (a>0) then SignLim:=1;
   if (a<-lim) and (a<0) then SignLim:=-1;
   if Abs(a)<lim then SignLim:=0;
end; //SignLim


function MaxReal(a,b:real):real;

begin
   if a>=b then
         MaxReal:=a
      else
         MaxReal:=b;
end; //MaxReal


function MinReal(a,b:real):real;

begin
   if a<=b then
         MinReal:=a
      else
         MinReal:=b;
end; //MinReal


function MaxInt(a,b:integer):integer;

begin
   if a>=b then
         MaxInt:=a
      else
         MaxInt:=b;
end; //MaxInt


function MinInt(a,b:integer):integer;

begin
   if a<=b then
         MinInt:=a
      else
         MinInt:=b;
end; //MinInt


procedure ChangeReal(var a,b:real);

var
   s:    real;

begin
   s:=a;
   a:=b;
   b:=s;
end; //ChangeReal


procedure ChangeInt(var a,b:integer);

var
   s:    integer;

begin
   s:=a;
   a:=b;
   b:=s;
end; //ChangeReal


function PowInt(a:real;b:integer):real;

var
   i:    integer;
   s:    real;

begin
   s:=1;
   for i:=1 to b do
      s:=s*a;
   PowInt:=s;
end; //PowInt


function PowReal(a,b:real):real;

begin
   PowReal:=exp(b*ln(a));
end; //PowReal


function TanD(alf:real):real;

begin
   alf:=alf/180*Pi;
   TanD:=Sin(alf)/Cos(alf);
end; //Tan


function ATan(dx,dy:real):real;

begin
   if dx>0 then
      if dy>=0 then
            ATan:=ArcTan(dy/dx)
         else
            ATan:=ArcTan(dy/dx)+2*Pi;
   if dx=0 then
      begin
         if dy>0 then ATan:=0.5*Pi;
         if dy=0 then ATan:=0;
         if dy<0 then ATan:=1.5*Pi;
      end;
   if dx<0 then ATan:=ArcTan(dy/dx)+Pi;

end; //ATan


function Interpol(x1,y1,x2,y2,x:real):real;

begin
   Interpol:=y1+(x-x1)/(x2-x1)*(y2-y1);
end; //Interpol



function DecConv(s:shortstring):shortstring;

var
   i,l:   integer;

begin
   l:=Length(s);
   for i:=1 to l do
      if s[i]=',' then s[i]:='.';
   DecConv:=s;
end; //DecConv



function DecConv2(s:shortstring):shortstring;

var
   i,l:       integer;
   s1:        shortstring;
   ch1,ch2:   char;

begin
   s1:=FloatToStr(1/2);
   if s1[2]='.' then
      begin
         ch1:='.';
         ch2:=',';
      end;
   if s1[2]=',' then
      begin
         ch1:=',';
         ch2:='.';
      end;
   l:=Length(s);
   for i:=1 to l do
      if s[i]=ch2 then s[i]:=ch1;
   DecConv2:=s;
end; //DecConv


function STrunc(s:shortstring;n:integer):shortstring;

var
   i: integer;
   t: shortstring;

begin
   t:=Copy(s,1,MinInt(n,Length(s)));
   //for i:=1 to MinInt(n,Length(s)) do
   //   t[i]:=s[i];
   for i:=Length(s)+1 to n do
      t:=t+' ';                                                    
   STrunc:=Copy(t,1,n);
end; //STrunc


procedure StringIn(s:string;var v:realnodetype);

var
   l,i,j1,j2,i2,n:                      integer;

begin
   s:=DecConv2(s);
   for i:=1 to maxnode do
      v[i]:=0;
   l:=Length(s);
   i:=1;
   n:=0;
   repeat
      begin
         repeat
            j1:=Pos(' ',s);
            j2:=Pos(Chr(9),s);
            i2:=Round(MinReal(j1,j2));
            if i2=0 then i2:=Round(MaxReal(j1,j2));
            if i2=1 then
               begin
                  Delete(s,1,1);
                  n:=n+1;
               end;
         until (i2<>1) or (n=l);
         if i2>1 then
            begin
               v[i]:=StrToFloat(Copy(s,1,i2-1));
               i:=i+1;
               Delete(s,1,i2);
               n:=n+i2;
            end;
         if i2=0 then v[i]:=StrToFloat(s);
      end;
   until (n=l) or (i2=0);
end; //StringIn


procedure StringIn2(s:string;var v:realnodetype);

var
   l,i,j1,j2,i2,n:                      integer;

begin
   s:=DecConv2(s);
   for i:=1 to maxnode do
      v[i]:=0;
   i:=1;
   repeat
      j1:=Pos('"',s);
      Delete(s,1,j1);
      j1:=Pos('"',s);
      v[i]:=StrToFloat(Copy(s,1,j1-1));
      Delete(s,1,j1);
      i:=i+1;
   until Pos(',',s)=0;;
end; //StringIn2


procedure StringIn3(s:string;var v:string10type);

var
   l,i,j1,j2,i2,n:                      integer;

begin
   s:=DecConv2(s);
   for i:=1 to 20 do
      v[i]:='';
   i:=1;
   repeat
      j1:=Pos('"',s);
      Delete(s,1,j1);
      j1:=Pos('"',s);
      v[i]:=Copy(s,1,j1-1);
      Delete(s,1,j1);
      i:=i+1;
   until Pos(',',s)=0;
end; //StringIn3



function StrBegin(s1,s2:shortstring):boolean;

//function checks whether the beginning of s1 is identical to s2

var
   b:   boolean;
   i:   integer;

begin
   b:=True;
   if Length(s1)<Length(s2)
      then b:=False
      else
         for i:=1 to Length(s2) do
            if s1[i]<>s2[i] then b:=False;
   StrBegin:=b;
end; //StrBegin



function CheckParallel(v1,v2:tvonal):boolean;

begin
   if Abs( 1 - Abs(((v1.x2-v1.x1)/(v1.y2-v1.y1)) / ((v2.x2-v2.x1)/(v2.y2-v2.y1))) ) < 0.01 then CheckParallel:=True
                                                                                           else CheckParallel:=False;
end; //CheckParallel



procedure LtoG(var dx,dy:real;alf:real);

var
      sx,sy:  real;

begin
   sx:=dx;
   sy:=dy;
   dx:=-sy*Sin(alf)+sx*Cos(alf);
   dy:=sx*Sin(alf)+sy*Cos(alf);
end; //LtoG




//*********************
//
//   Input Routines
//
//*********************



procedure ReadInpData(var nbeam,nnode,nsup:integer;var nodes:nodeproptype;
                      var beamtop:beamtopoltype;var supports:nodeproptype;
                      var supploads:reallcnodeproptype;
                      var E,alft:real;var Av,Iv:realbeamtype;var nlc:integer;
                      var nloads:intlcbeamtype;var loadtyps:intlcbeamloadtype;
                      var loadprops:multiloadproptype;var tonew,toold:intnodetype);

var
   i,j,ilc,ib,jb,iorig,jorig,nprof:  integer;
   x1,x2,y1,y2,t,A,I1,I2,cg:  real;
   f:  text;
   datadir,name,nam,nam2,s:  string;
   v:  realnodetype;
   nonstruct,braceconnect,ch: boolean;


begin
   //nodes
   nnode:=Length(mainf.geom.aStatCsp);
   for i:=1 to nnode do
      begin
         nodes[i,1]:=mainf.geom.aStatCsp[i-1].x;
         nodes[i,2]:=-mainf.geom.aStatCsp[i-1].y;  //minus sign is to consider the difference between the coord systems
      end;
   nbeam:=Length(mainf.geom.aStatVaz);

   //renumbering nodes to reduce band-width
   for i:=1 to nnode do
      tonew[i]:=i;
//{
   repeat
      ch:=True;
      for i:=1 to nnode-1 do
         begin
            if nodes[i,1]>nodes[i+1,1] then
               begin
                  ChangeReal(nodes[i+1,1],nodes[i,1]);
                  ChangeReal(nodes[i+1,2],nodes[i,2]);
                  //ChangeReal(nodes[i+1,3],nodes[i,3]);
                  //ChangeReal(nodes[i+1,4],nodes[i,4]);
                  ChangeInt(tonew[i+1],tonew[i]);
                  ch:=False;
               end;
         end;
   until ch;                   //  }
   for i:=1 to nnode do
      toold[tonew[i]]:=i;


   //beam elements
   for i:=1 to nbeam do
      begin
         x1:=mainf.geom.aStatVaz[i-1].v.x1;
         y1:=-mainf.geom.aStatVaz[i-1].v.y1;  //minus sign due to cord syst transform
         j:=0;
         repeat
            j:=j+1;
         until (Abs(nodes[j,1]-x1)<0.000001) and (Abs(nodes[j,2]-y1)<0.000001);
         beamtop[i,1]:=j;
         x2:=mainf.geom.aStatVaz[i-1].v.x2;
         y2:=-mainf.geom.aStatVaz[i-1].v.y2;
         j:=0;
         repeat
            j:=j+1;
         until (Abs(nodes[j,1]-x2)<0.000001) and (Abs(nodes[j,2]-y2)<0.000001);
         beamtop[i,2]:=j;
      end;

   //supports
   nsup:=0;
   for i:=1 to nnode do
      begin
         if mainf.geom.aStatCsp[i-1].CsuklosTam then
            begin
               nsup:=nsup+1;
               supports[nsup,1]:=toold[i];
               supports[nsup,2]:=FixSpring;
               supports[nsup,3]:=FixSpring;
               supports[nsup,4]:=0;
            end;
         if mainf.geom.aStatCsp[i-1].GorgosTam then
            begin
               nsup:=nsup+1;
               supports[nsup,1]:=toold[i];
               supports[nsup,2]:=0;
               supports[nsup,3]:=FixSpring;
               supports[nsup,4]:=0;
            end;
      end;

   //material props
   E:=210000;    //modulus of elasticity
   alft:=1.2e-5;   //coefficient of thermal expansion


   //cross-section properties

   for ib:=1 to nbeam do
      begin
         iorig:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
         name:=mainf.geom.aRud[iorig-1].profile.profile;

         //section props for one profile
          nprof:=0;
          repeat
             nprof:=nprof+1;
          until aProfile[nprof-1].profile=name;
          A:=aProfile[nprof-1].Ag;
          I1:=aProfile[nprof-1].Iyg;
          I2:=aProfile[nprof-1].Izg;
          cg:=aProfile[nprof-1].yCg;


          //checking if the element is a structural element or not
          nonstruct:=False;
          if mainf.geom.aStatVaz[ib-1].Stagbar then nonstruct:=True;
          if mainf.geom.aStatVaz[ib-1].PurlinSupport and mainf.chkFugg.Checked then nonstruct:=True;
          if mainf.geom.aStatVaz[ib-1].PurlinSupport and mainf.chkMer.Checked then
             begin
                braceconnect:=False;
                for jb:=1 to nbeam do
                   if jb<>ib then
                      begin
                         if (beamtop[jb,1]=beamtop[ib,1]) or (beamtop[jb,1]=beamtop[ib,2])
                            or (beamtop[jb,2]=beamtop[ib,1]) or (beamtop[jb,2]=beamtop[ib,2]) then
                               begin
                                  jorig:=mainf.geom.aStatvaz[jb-1].RudIdx+1;
                                  if mainf.geom.ARud[jorig-1].tipus in [rBelsoOszlop] then
                                     braceconnect:=True;
                               end;
                      end;
                if not braceconnect then nonstruct:=True;
             end;

         //defining final properties
         if nonstruct
            then
               begin
                  Av[ib]:=RigidA;
                  Iv[ib]:=RigidI;
               end
            else
               begin

         if mainf.geom.ARud[iorig-1].tipus in [rAlsoOv,rFelsoOv] then
            begin
               if mainf.geom.ARud[iorig-1]._U then
                  begin
                     Av[ib]:=1.8*2*A;
                     Iv[ib]:=1.8*2*I1;
                  end
               else
                  begin
                     Av[ib]:=2*A;
                     Iv[ib]:=2*I1;
                  end;
            end;
         if mainf.geom.ARud[iorig-1].tipus in [rBak,rKulsoOszlop] then
            begin
               Av[ib]:=A;
               Iv[ib]:=I2;
            end;
         if mainf.geom.ARud[iorig-1].tipus in [rBelsoOszlop,rRacs] then
            begin
               if mainf.geom.ARud[iorig-1].dupla then
                  begin
                     Av[ib]:=2*A;
                     Iv[ib]:=2*(I2+A*cg*cg);   //'cg' a s?lypont t?vols?ga a gerinct?l!!!
                  end
               else
                  begin
                     Av[ib]:=A;
                     Iv[ib]:=I2;
                  end
            end;

               end;
      end; //for ib



   //loads
   nlc:=4;  //nr of load cases
   for ilc:=1 to maxcase do
      for ib:=1 to maxbeam do
         nloads[ilc,ib]:=0;
   for ilc:=1 to maxcase do
      for ib:=1 to maxbeam do
         for j:=1 to maxload do
            loadtyps[ilc,ib,j]:=0;
   for ilc:=1 to maxcase do
      for ib:=1 to maxnode do
         for j:=1 to 4 do
            supploads[ilc,ib,j]:=0;

       //tesztel?sre valami
       //1. loadcase
       //a default szerkezet eset?n az als? ?v?n egy koncentr?lt er?
   i:=30;  //az i.-ik r?don van teher
   nloads[1,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[1,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[1,i,1,1]:=0;//50000;  //az er? ?rt?ke [N]
   loadprops[1,i,1,3]:=0.3;   //az er? relat?v helyzete a r?don bel?l (0..1)
   //supploads[1,1,1]:=5; //5 mm v?zszintesen jobbra
   //supploads[1,2,2]:=5; //5 mm f?gg?legesen lefel?

       //2. loadcase
       //a default szerkezet eset?n az als? ?v?n egy koncentr?lt er?
   i:=29;  //az i.-ik r?don van teher
   nloads[2,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[2,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[2,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[2,i,1,3]:=0.0;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=39;  //az i.-ik r?don van teher
   nloads[2,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[2,i,1]:=5; //a teher t?pusa: koncentr?lt nyomat?k
   loadprops[2,i,1,1]:=50e6;  //a nyomat?k ?rt?ke [Nmm]
   loadprops[2,i,1,3]:=0.9;   //a nyomat?k relat?v helyzete a r?don bel?l (0..1)

    //3. loadcase
   //a default szerkezet eset?n a fels? ?v?n mindenhol egy koncentr?lt er?
   i:=14;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=21;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=28;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=33;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)

   i:=44;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=51;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=58;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=64;  //az i.-ik r?don van teher
   nloads[3,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[3,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[3,i,1,1]:=50000;  //az er? ?rt?ke [N]
   loadprops[3,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)

    //4. loadcase
   //a default szerkezet eset?n a fels? ?v?n mindenhol egy koncentr?lt er?
   i:=14;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=21;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=28;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=33;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)

   i:=44;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=51;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=58;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)
   i:=64;  //az i.-ik r?don van teher
   nloads[4,i]:=1;  //az i.-ik r?don a tehereset komponenseinek darabsz?ma
   loadtyps[4,i,1]:=1; //a teher t?pusa: koncentr?lt er? (mer?leges a r?dra)
   loadprops[4,i,1,1]:=5000;  //az er? ?rt?ke [N]
   loadprops[4,i,1,3]:=0.5;   //az er? relat?v helyzete a r?don bel?l (0..1)



//ideiglenesen: az geom adatok ki?rat?sa
   Assign(f,'geom.txt');
   ReWrite(f);
   Writeln(f,'CSOM?PONTI KOORDIN?T?K');
   Writeln(f,'sorsz?m, x [mm],   y[mm]');
   for i:=1 to nnode do
      Writeln(f,i:5,nodes[i,1]:12:5,nodes[i,2]:12:5);
   Writeln(f);   
   Writeln(f,'R?DELEMEK TOPOL?GI?JA, KERESZTMETSZETI ADATAI');
   Writeln(f,'sorsz, kezd?p, v?gp,  ter?let,   inercia,  rug.mod.');
   Writeln(f,'                      [mm2]       [mm4]     [N/mm2]'); 
   for i:=1 to nbeam do
      Writeln(f,i:5,beamtop[i,1]:5,beamtop[i,2]:5,Av[i]:12:1,Iv[i]:12:0,E:12:0);
   Writeln(f);
   Writeln(f,'T?MASZOK (rug?kkal defini?lva)');
   Writeln(f,'sorsz, csom?pont,       r_x,           r_y,          r_fi');
   Writeln(f,'                       [N/mm]         [N/mm]         [N/mm]');
   for i:=1 to nsup do
      Writeln(f,i:5,round(supports[i,1]):10,supports[i,2]:15:0,supports[i,3]:15:0,supports[i,4]:15:0);
   Close(f);


end; //ReadInpData






//**************************************************
//
//Subroutines for Solving Linear System of Equations
//
//**************************************************



procedure SeparCholesky(n,w:integer;A:stifmattype;var L:stifmattype);

//to make Cholesky separation

//n        nro of rows
//w        bandwidth
//A        matrix to be separated
//L        lower triangle matrix

var
   i,j,k:        integer;
   s:            double;

begin
   for i:=1 to w-1 do
      for j:=1 to (w-i) do
         L[i,j]:=0;

   for i:=1 to n do
      for j:=(i+1-w) to i do
         if j>0 then
            begin
               if i>j then
                  begin
                     s:=0;
                     for k:=(i+1-w) to (j-1) do
                        s:=s+L[i,k-i+w]*L[j,k-j+w];
                     L[i,j-i+w]:=(A[i,j-i+w]-s)/L[j,w];
                  end;
               if i=j then
                  begin
                     s:=0;
                     for k:=(i+1-w) to (i-1) do
                        s:=s+Sqr(L[i,k-i+w]);
                     if (A[i,w]-s)>1e-10 then L[i,w]:=Sqrt(A[i,w]-s)
                                         else exit;
                  end;
            end;
end; //SeparCholesky


procedure BackSubstitute(n,w,nlc:integer;L:stifmattype;b:reallcdoftype;
                         var x:reallcdoftype);

//to calculate x vector from L

//n        nro of rows
//w        bandwidth
//nlc      nro of load cases
//L        lower triangle matrix
//b        matrix on the (multiple) right-hand side
//x        the vector of vector of unknowns

var
   ic,i,k:   integer;
   s:        real;

begin
   for ic:=1 to nlc do
      begin

         for i:=1 to n do
            begin
               s:=0;
               for k:=(i+1-w) to (i-1) do
                  if k>0 then
                     s:=s+L[i,k-i+w]*x[ic,k];
               x[ic,i]:=(b[ic,i]-s)/L[i,w];
            end;

         for i:=n downto 1 do
            begin
               s:=0;
               for k:=(i+1) to n do
                  if (k-i)<w then
                     s:=s+L[k,i-k+w]*x[ic,k];
               x[ic,i]:=(x[ic,i]-s)/L[i,w]
            end;

      end; //for ic
end; //BackSubstitute


procedure EquaSolve(n,w,nlc:integer;A:stifmattype;b:reallcdoftype;var x:reallcdoftype);

var
   L:     stifmattype;

begin
   SeparCholesky(n,w,A,L);
   BackSubstitute(n,w,nlc,L,b,x);
end; //EquaSolve



//******************************
//
//Routins of Displacement Method
//
//******************************


function BandCalc(nb:integer;beamtop:beamtopoltype):integer;

//to calculate bandwidth of the global stiffness matrix

var
   i,nw:  integer;

begin
   nw:=0;
   for i:=1 to nb do
      nw:=MaxInt(Abs(beamtop[i,1]-beamtop[i,2]),nw);
   BandCalc:=(nw+1)*3;
end; //BandCalc


procedure BeamProps(nn,nb:integer;nodes:nodeproptype;beamtop:beamtopoltype;
                    var Lv,alfv:realbeamtype);

//to calculate cross-sectional properties and angles for the beam elements

var
   i:  integer;
   dx,dy:  real;

begin
   for i:=1 to nb do
      begin
         dx:=nodes[beamtop[i,2],1]-nodes[beamtop[i,1],1];
         dy:=nodes[beamtop[i,2],2]-nodes[beamtop[i,1],2];
         Lv[i]:=Sqrt(dx*dx+dy*dy);
         alfv[i]:=ATan(dx,dy);
         //alfv[i]:=ATan(dx,dy)*180/Pi;
      end;
end; //BeamProps


procedure CalcKlocal(E,A,I,L:real;var Kloc:stifloctype);

//to calculate local stiffness matrix for a beam element

var
   k,j:           integer;
   k1,k2,k3,k4:   real;

begin
   k1:=E*A/L;
   k2:=2*E*I/L;   // 2EI/L
   k3:=k2*3/L;    // 6EI/L^2
   k4:=k3*2/L;   // 12EI/L^3
   for k:=1 to 6 do
      for j:=k to 6 do
         Kloc[k,j]:=0;
   Kloc[1,1]:=k1;
   Kloc[1,4]:=-k1;
   Kloc[2,2]:=k4;
   Kloc[2,3]:=k3;
   Kloc[2,5]:=-k4;
   Kloc[2,6]:=k3;
   Kloc[3,3]:=2*k2;
   Kloc[3,5]:=-k3;
   Kloc[3,6]:=k2;
   Kloc[4,4]:=k1;
   Kloc[5,5]:=k4;
   Kloc[5,6]:=-k3;
   Kloc[6,6]:=2*k2;
   for k:=1 to 6 do
      for j:=1 to (k-1) do
         Kloc[k,j]:=Kloc[j,k];
end;  //CalcKlocal


procedure Transf(alf:real;var T:stifloctype);

//to calculate transformation matrix for coordinate rotation

var
   i,j:  integer;
   c,s:  real;

begin
   for i:=1 to 6 do
      for j:=1 to 6 do
         T[i,j]:=0;
   c:=cos(alf);
   s:=sin(alf);
   T[1,1]:=c;
   T[1,2]:=-s;
   T[2,1]:=s;
   T[2,2]:=c;
   T[3,3]:=1;
   T[4,4]:=c;
   T[4,5]:=-s;
   T[5,4]:=s;
   T[5,5]:=c;
   T[6,6]:=1;
end; //Transf


procedure Ktransf(T:stifloctype;var K:stifloctype);


//to calculate local stiffness matrix transformed to the global coordinate system

var
   i,j,m: integer;
   s:     real;
   TK:    stifloctype;

begin
   for i:=1 to 6 do
      for j:=1 to 6 do
         begin
            s:=0;
            for m:=1 to 6 do
               s:=s+T[i,m]*K[m,j];
            TK[i,j]:=s;
         end;
   for i:=1 to 6 do
      for j:=1 to 6 do
         begin
            s:=0;
            for m:=1 to 6 do
               s:=s+TK[i,m]*T[j,m];
            K[i,j]:=s;
         end;
end; //Ktransf


procedure AddSupports(nsup:integer;supports:nodeproptype;var Kglob:doubledofdoftype);

//to add springs to global stiffness matrix to conside supports

//nsup       number of supported nodes
//supports   vector (matrix) of support data
//               number of rows = nsup
//               in a row: node nr, x_dir, y_dir, rotation (spring stiffness)

var
   i,j,n:   integer;

begin
   for i:=1 to nsup do
      begin
         n:=(Round(supports[i,1])-1)*3;
         for j:=1 to 3 do
            Kglob[n+j,n+j]:=Kglob[n+j,n+j]+supports[i,j+1];
      end;
end; //Addsupports


procedure Kglobal(nb,nn,nw,nsup:integer;beamtop:beamtopoltype;supports:nodeproptype;
                  Av,Iv,Lv,alfv:realbeamtype;E:real;
                  var Kglob:stifmattype);

//to add beam stiffnesses to the stiffness matrix

//nb      nro of beam elements
//nn      nro of nodes
//nw      bandwidth (actually: half bandwidth) of stiffness matrix
//beamtop vector of beam topology [start node nr, end node nr]
//Iv      vector of inertias
//Lv      vector of beam lengths
//alfv    vector of beam angles
//E       modulus of elasticity
//Kglob   global stiffness matrix (bandmatrix)
//Ks      same as Kglob, but it is a square matrix
//            (later, Ks can be eliminated if necessary)

var
   i,j,k,n1,n2:  integer;
   Kloc,T:       stifloctype;
   Ks:           doubledofdoftype;

begin
   //zero out Ks, Kglob
   for i:=1 to maxdof do
      for j:=1 to maxdof do
         Ks[i,j]:=0;
   for i:=1 to maxdof do
      for j:=1 to maxband do
         Kglob[i,j]:=0;

   //compile Ks
   for i:=1 to nb do
      begin
         CalcKlocal(E,Av[i],Iv[i],Lv[i],Kloc);
         Transf(alfv[i],T);
         Ktransf(T,Kloc);
         n1:=beamtop[i,1];
         n2:=beamtop[i,2];
         n1:=(n1-1)*3;
         n2:=(n2-1)*3;
         for j:=1 to 3 do
            for k:=1 to 3 do
               begin
                  Ks[n1+j,n1+k]:=Ks[n1+j,n1+k]+Kloc[j,k];
                  Ks[n1+j,n2+k]:=Ks[n1+j,n2+k]+Kloc[j,k+3];
                  Ks[n2+j,n1+k]:=Ks[n2+j,n1+k]+Kloc[j+3,k];
                  Ks[n2+j,n2+k]:=Ks[n2+j,n2+k]+Kloc[j+3,k+3];
               end;
      end;

      //add supports
      AddSupports(nsup,supports,Ks);


      //transform Ks to bandmatrix Kglob
      for j:=1 to (3*nn) do
         for k:=1 to nw do
            if j>(nw-k) then Kglob[j,k]:=Ks[j,j-(nw-k)]
                        else Kglob[j,k]:=0;
end; //Kglobal



procedure Calc_qloc(L,Area,E,alft:real;qltyp:integer;vload:vloadtype;var qloc:real6type);

//to calculate local load vector for a beam element

//Feltev?sek:
//  vload[1] koncentr?lt teher intenzit?sa
//           vagy megoszl? teher kezd?intenzit?sa (kezd?ponthoz k?zelebbi intenzit?s)
//  vload[2] (koncentr?lt teher eset?n inakt?v adat)
//           megoszl? teher v?gintenzit?sa (kezd?pontt?l t?volabbi intenzit?s)
//  vload[3] a kezd?ponthoz k?pesti relat?v teherhelyzet (koncentr?lt er?/nyom. eset?n)
//           vagy a kezd?ponthoz k?zelebbi teherv?g relat?v helyzete a kezd?ponthoz k?pest (megoszl? teher eset?n)
//  vload[4] (koncentr?lt er?/nyom. eset?n inakt?v adat)
//           a kezd?pontt?l t?volabbi teherv?g relat?v helyzete a kezd?opnthoz k?pest (megoszl? teher eset?n)


var
   a,b,c,s,F,M,q1,q2,m1,m2,dt:   real;

begin
   case qltyp of
      1: begin  //conc force perp
            F:=vload[1];
            a:=vload[3]*L;
            c:=L-a;
            qloc[1]:=0;
            qloc[2]:=F*PowInt(c,2)*(3*a+c)/PowInt(L,3);
            qloc[3]:=F*a*PowInt(c/L,2);
            qloc[4]:=0;
            qloc[5]:=F*PowInt(a,2)*(3*c+a)/PowInt(L,3);
            qloc[6]:=-F*c*PowInt(a/L,2);
         end;
      2: begin  //conc force parallel
            F:=vload[1];
            a:=vload[3]*L;
            c:=L-a;
            qloc[1]:=F*c/L;
            qloc[2]:=0;
            qloc[3]:=0;
            qloc[4]:=F*a/L;
            qloc[5]:=0;
            qloc[6]:=0;
         end;
      3: begin  //dist load perp
            q1:=vload[1];
            q2:=vload[2];
            a:=vload[3]*L;
            c:=(1-vload[4])*L;
            b:=L-a-c;
            qloc[1]:=0;
            s:=PowInt(b,3)*(7*q1+3*q2) + 10*PowInt(c,3)*(q1+q2) + 5*a*PowInt(b,2)*(3*q1+q2) ;
            s:=s + 5*c*PowInt(b,2)*(5*q1+3*q2) + 30*a*PowInt(c,2)*(q1+q2) ;
            s:=s + 30*b*PowInt(c,2)*(q1+q2) + 20*a*b*c*(2*q1+q2);
            qloc[2]:=b/20/PowInt(L,3)*s;
            s:=PowInt(b,3)*(3*q1+2*q2) + 5*a*PowInt(b,2)*(3*q1+q2) + 10*c*PowInt(b,2)*(q1+q2) ;
            s:=s + 30*a*PowInt(c,2)*(q1+q2) + 10*b*PowInt(c,2)*(q1+2*q2) + 20*a*b*c*(2*q1+q2);
            qloc[3]:=b/60/PowInt(L,2)*s;
            qloc[4]:=0;
            s:=PowInt(b,3)*(3*q1+7*q2) + 10*PowInt(a,3)*(q1+q2) + 5*a*PowInt(b,2)*(3*q1+5*q2) ;
            s:=s + 5*c*PowInt(b,2)*(q1+3*q2) + 30*c*PowInt(a,2)*(q1+q2) ;
            s:=s + 30*b*PowInt(a,2)*(q1+q2) + 20*a*b*c*(q1+2*q2);
            qloc[5]:=b/20/PowInt(L,3)*s;
            s:=PowInt(b,3)*(2*q1+3*q2) + 10*a*PowInt(b,2)*(q1+q2) +  5*c*PowInt(b,2)*(q1+3*q2) ;
            s:=s + 30*c*PowInt(a,2)*(q1+q2) + 10*b*PowInt(a,2)*(2*q1+q2) + 20*a*b*c*(q1+2*q2);
            qloc[6]:=-b/60/PowInt(L,2)*s;
         end;
      4: begin  //dist load parallel
            q1:=vload[1];
            q2:=vload[2];
            a:=vload[3]*L;
            c:=(1-vload[4])*L;
            b:=L-a-c;
            s:=0.5*(q1+q2)*b;
            //qloc[1]:=s*(0.5*b+c)/L;
            qloc[1]:=1/6*b*(q2*b+2*q1*b+3*c*q1+3*c*q2)/L;
            qloc[2]:=0;
            qloc[3]:=0;
            //qloc[4]:=s*(0.5*b+a)/L;
            qloc[4]:=1/6*b*(2*q2*b+q1*b+3*a*q1+3*a*q2)/L;
            qloc[5]:=0;
            qloc[6]:=0;
         end;
      5: begin  //conc moment
            M:=vload[1];
            a:=vload[3]*L;
            c:=L-a;
            qloc[1]:=0;
            qloc[2]:=-M*6*a*c/PowInt(L,3);
            qloc[3]:=-M*c*(2*a-c)/PowInt(L,2);
            qloc[4]:=0;
            qloc[5]:=M*6*a*c/PowInt(L,3);
            qloc[6]:=-M*a*(2*c-a)/PowInt(L,2);
         end;
      6: begin  //dist moment
            m1:=vload[1];
            m2:=vload[2];
            a:=vload[3]*L;
            c:=(1-vload[4])*L;
            b:=L-a-c;
            qloc[1]:=0;
            qloc[2]:=-( (b*b+6*a*c)*(m1+m2) + 2*a*b*(2*m1+m2) + 2*b*c*(m1+2*m2) ) *b / (2*PowInt(L,3));
            qloc[3]:=-( (12*a*c-6*c*c)*(m1+m2) + (b*b+4*b*c)*(-m1+m2) + 4*a*b*(2*m1+m2) ) *b / (12*L*L);
            qloc[4]:=0;
            qloc[5]:=-qloc[2];
            qloc[6]:=-( (12*a*c-6*a*a)*(m1+m2) + (b*b+4*b*a)*(m1-m2) + 4*c*b*(m1+2*m2) ) *b / (12*L*L);
         end;
      7: begin  //temperature change
            dt:=vload[1];
            qloc[1]:=-alft*dt*E*Area;
            qloc[2]:=0;
            qloc[3]:=0;
            qloc[4]:=alft*dt*E*Area;
            qloc[5]:=0;
            qloc[6]:=0;
         end;
      8: begin  //support displacement
            qloc[1]:=0;
            qloc[2]:=0;
            qloc[3]:=0;
            qloc[4]:=0;
            qloc[5]:=0;
            qloc[6]:=0;
         end;
   end; //case
end; //Calc_qloc



procedure Transf_q(T:stifloctype;var qloc:real6type);

//to calculate local load vector transformed to the global coordinate system

var
   i,j:   integer;
   s:     real;
   Tq:    real6type;

begin
   for i:=1 to 6 do
      begin
         s:=0;
         for j:=1 to 6 do
            s:=s+T[i,j]*qloc[j];
         Tq[i]:=s;
      end;
   for i:=1 to 6 do
      qloc[i]:=Tq[i];

end; //Transf_q


procedure Calcqglob(nlc,nsup:integer;nloads:intlcbeamtype;loadtyps:intlcbeamloadtype;
                    loadprops:multiloadproptype;nb:integer;beamtop:beamtopoltype;
                    Lv,Av,alfv:realbeamtype;supports:nodeproptype;
                    supploads:reallcnodeproptype;E,alft:real;
                    var qglob:reallcdoftype);

//to calculate global load vector

//nlc         nr of load cases
//nsup        nr of supports
//nloads      1st dim: 1..(nr of load cases)
//            2nd dim: 1..(nr of beam elements)
//            value: nr of sub-loads
//loadtyps    1st dim: 1..(nr of load cases)
//            2nd dim: 1..(nr of beam elements)
//            3rd dim: 1..(nr of sub-loads)
//            value: type of subload
//loadprops   1st dim: 1..(nr of load cases)
//            2nd dim: 1..(nr of beam elements)
//            3rd dim: 1..(nr of sub-loads)
//            value: properties of subloads (which is a vector of vloadtype)
//nb          beam element identifier
//beamtop     beam topology (start node, end node)
//Lv,Av       length and cs area of beam elements
//supports    1st dim: 1..(nr of nodes)
//            2nd dim: 1..4 (1:node nr, 2:x spring, 3:y spring, 4:rotational spring)
//supploads   1st dim: 1..(nr of load cases)
//            2nd dim: 1..(nr of nodes)
//            3rd dim: 1..3 (1:x displ, 2:y displ, 3:rotational displ)
//alfv        angle of beam element (wrt positive global x axis)
//E,alft      Young's modulus, coefficient of linear thermal expansion
//qglob       1st dim: 1..(nr of load cases)
//            2nd dim: 1..(nr of nodal dofs)
//            value: the value of the given nodal force/moment



var
   i,j,k,ic,nl,n1,n2:  integer;
   ql,qloc:            real6type;
   T:                  stifloctype;

begin
   //zero out qglob
   for ic:=1 to nlc do
      for i:=1 to maxdof do
         qglob[ic,i]:=0;

   //compile qglob
   for ic:=1 to nlc do
      begin

         for i:=1 to nb do
            begin
               for j:=1 to 6 do
                  qloc[j]:=0;
               nl:=nloads[ic,i];
               for j:=1 to nl do
                  begin
                     Calc_qloc(Lv[i],Av[i],E,alft,loadtyps[ic,i,j],loadprops[ic,i,j],ql);
                     for k:=1 to 6 do
                        qloc[k]:=qloc[k]+ql[k];
                  end;
               if nl>0 then
                  begin
                     Transf(alfv[i],T);
                     Transf_q(T,qloc);
                     n1:=beamtop[i,1];
                     n2:=beamtop[i,2];
                     n1:=(n1-1)*3;
                     n2:=(n2-1)*3;
                     for j:=1 to 3 do
                        begin
                           qglob[ic,n1+j]:=qglob[ic,n1+j]+qloc[j];
                           qglob[ic,n2+j]:=qglob[ic,n2+j]+qloc[j+3];
                        end;
                  end;
            end;  //for i (beam elements)

         //add loads from support displacements
         for i:=1 to nsup do
            begin
               n1:=Round(supports[i,1]);
               qglob[ic,3*(n1-1)+1]:=qglob[ic,3*(n1-1)+1]+supploads[ic,i,1]*supports[i,2];
               qglob[ic,3*(n1-1)+2]:=qglob[ic,3*(n1-1)+2]+supploads[ic,i,2]*supports[i,3];
               qglob[ic,3*(n1-1)+3]:=qglob[ic,3*(n1-1)+3]+supploads[ic,i,3]*supports[i,4];
            end; //for i (supports)

      end; //for ic (load cases)


end; //Calcqglob


procedure CalcIntForces(nloads:intlcbeamtype;loadtyps:intlcbeamloadtype;loadprops:multiloadproptype;
                        nlc,nb:integer;beamtop:beamtopoltype;Lv,Av,Iv,alfv:realbeamtype;E,alft:real;
                        xglob:reallcdoftype;npoints:integer;xs:realpointstype;
                        var Nel,Tel,Mel:reallcpointstype);

//to calculate moment (stress resultant) in a given cross-section of a beam

var
   i,j,k,n1,n2,nl,ip,ilc:           integer;
   s,L,a,b,c,F,M,q1,q2,m1,m2,dt:    real;
   x,nx,tx,mx,R1,R2:                real;
   intf:                            real3type;
   xloc,Txloc,fx,fq,fs,fend:        real6type;
   T,Kloc:                          stifloctype;

begin
   //preparation
   L:=Lv[nb];
   Transf(alfv[nb],T);
   CalcKlocal(E,Av[nb],Iv[nb],Lv[nb],Kloc);

   //start dof, end dof
   n1:=(beamtop[nb,1]-1)*3;
   n2:=(beamtop[nb,2]-1)*3;

   for ilc:=1 to nlc do
      begin

         nl:=nloads[ilc,nb];

         //local displacement vector in global coordinate system
         for i:=1 to 3 do
            begin
               xloc[i]:=xglob[ilc,n1+i];
               xloc[i+3]:=xglob[ilc,n2+i];
            end;

         //transform displacements from global to local coordinate system
         for i:=1 to 6 do
            begin
               s:=0;
               for j:=1 to 6 do
                  s:=s+T[j,i]*xloc[j];    //T_transp ? xloc
               Txloc[i]:=s;
            end;
         for i:=1 to 6 do
            xloc[i]:=Txloc[i];

         //end moment/forces from nodal displacments
         for i:=1 to 6 do
            begin
               s:=0;
               for j:=1 to 6 do
                  s:=s+Kloc[i,j]*xloc[j];    //Kloc ? xloc
               fx[i]:=s;
            end;

         //end forces/moments from loads
         for k:=1 to 6 do
            fq[k]:=0;
         for j:=1 to nl do
            begin
               Calc_qloc(Lv[nb],Av[nb],E,alft,loadtyps[ilc,nb,j],loadprops[ilc,nb,j],fs);
               for k:=1 to 6 do
                  fq[k]:=fq[k]-fs[k];
            end;

         //end forces/moments
         for i:=1 to 6 do
            fend[i]:=fx[i]+fq[i];

         //cycle for the points begins
         for ip:=1 to npoints do
            begin
                x:=xs[ip];

               //internal forces/moments
               for i:=1 to 3 do
                  intf[i]:=0;
               x:=x*L;
               for j:=1 to nl do
                  begin
                     case loadtyps[ilc,nb,j] of
                        1: begin  //conc force perp
                              F:=loadprops[ilc,nb,j,1];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=L-a;
                              nx:=0;
                              if x<=a then tx:=0;
                              if x>a then tx:=F;
                              if x<=a then mx:=F*c/L*x;
                              if x>a then mx:=F*c/L*x-F*(x-a);
                           end;
                        2: begin  //conc force parallel
                              F:=loadprops[ilc,nb,j,1];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=L-a;
                              if x<=a then nx:=0;
                              if x>a then nx:=F;
                              tx:=0;
                              mx:=0;
                           end;
                        3: begin  //dist load perp
                              q1:=loadprops[ilc,nb,j,1];
                              q2:=loadprops[ilc,nb,j,2];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=(1-loadprops[ilc,nb,j,4])*L;
                              b:=L-a-c;
                              s:=0.5*(q1+q2)*b;
                              R1:=(q2*b/2*(b/3+c)+q1*b/2*(2*b/3+c))/L;
                              R2:=(q2*b/2*(a+2/3*b)+q1*b/2*(a+b/3))/L;
                              R2:=s-R1;
                              nx:=0;
                              if x<=a then tx:=0;
                              if (x>a) and (x<(a+b)) then tx:=0.5*(q1+( q1+(q2-q1)*(x-a)/b ))*(x-a);
                              if x>=(a+b) then tx:=s;
                              if x<=a then mx:=R1*x;
                              if (x>a) and (x<(a+b)) then mx:=x*R1-q2/6/b*PowInt(x-a,3)-q1/2*PowInt(x-a,2)+q1/6/b*PowInt(x-a,3);
                              if x>=(a+b) then mx:=R2*(L-x);
                           end;
                        4: begin  //dist load parallel
                              q1:=loadprops[ilc,nb,j,1];
                              q2:=loadprops[ilc,nb,j,2];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=(1-loadprops[ilc,nb,j,4])*L;
                              b:=L-a-c;
                              s:=0.5*(q1+q2)*b;
                              if x<=a then nx:=0;
                              if (x>a) and (x<(a+b)) then nx:=0.5*(q1+( q1+(q2-q1)*(x-a)/b ))*(x-a);
                              if x>=(a+b) then nx:=s;
                              tx:=0;
                              mx:=0;
                           end;
                        5: begin  //conc moment
                              M:=loadprops[ilc,nb,j,1];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=L-a;
                              nx:=0;
                              tx:=0;
                              if x<=a then mx:=-M/L*x;
                              if x>a then mx:=-M/L*x+M;
                              end;
                        6: begin  //dist moment
                              m1:=loadprops[ilc,nb,j,1];
                              m2:=loadprops[ilc,nb,j,2];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=(1-loadprops[ilc,nb,j,4])*L;
                              b:=L-a-c;
                              s:=0.5*(m1+m2)*b;
                              nx:=0;
                              tx:=0;
                              if x<=a then mx:=-s/L*x;
                              if (x>a) and (x<(a+b)) then mx:=-s/L*x+0.5*(m1+( m1+(m2-m1)*(x-a)/b ))*(x-a);
                              if x>=(a+b) then mx:=-s/L*x+s;
                           end;
                        7: begin  //temperature change
                              dt:=loadprops[ilc,nb,j,1];
                              nx:=0;
                              tx:=0;
                              mx:=0;
                           end;
                        8: begin  //support displacement
                              nx:=0;
                              tx:=0;
                              mx:=0;
                           end;
                     end; //case
                     intf[1]:=intf[1]+nx;
                     intf[2]:=intf[2]+tx;
                     intf[3]:=intf[3]+mx;
                  end; //for j

               //store N-T-M
                 //note, minus sign at N and T is to cosider engineering sign convention
               Nel[ilc,ip]:=-(fend[1]+intf[1]);
               Tel[ilc,ip]:=-(fend[2]+intf[2]);
               Mel[ilc,ip]:=(fend[3]*(L-x)/L-fend[6]*x/L+intf[3]);

            end; //for ip

      end; //for ilc

end; //CalcIntForces



procedure CalcReactions(nloads:intlcbeamtype;loadtyps:intlcbeamloadtype;loadprops:multiloadproptype;
                        nlc,nsup:integer;supports:nodeproptype;supploads:reallcnodeproptype;
                        xglob:reallcdoftype;
                        var reac:reallcnodeproptype);

//to calculate reaction forces of the structure

var
   ilc,isup,nod,j:           integer;

begin
   //zero out
   for ilc:=1 to maxcase do
      for isup:=1 to maxnode do
         for j:=1 to 4 do
            reac[ilc,isup,j]:=0;

   //calc reactions
   for ilc:=1 to nlc do
      begin
         for isup:=1 to nsup do
            begin
               nod:=Round(supports[isup,1]);
               reac[ilc,isup,4]:=nod;
               reac[ilc,isup,1]:=supports[isup,2]*(xglob[ilc,(nod-1)*3+1]-supploads[ilc,isup,1]);
               reac[ilc,isup,2]:=supports[isup,3]*(xglob[ilc,(nod-1)*3+2]-supploads[ilc,isup,2]);
               reac[ilc,isup,3]:=supports[isup,4]*(xglob[ilc,(nod-1)*3+3]-supploads[ilc,isup,3]);
            end; //isup

      end; //for ilc

end; //CalcReactions


procedure CalcDeformations(nloads:intlcbeamtype;loadtyps:intlcbeamloadtype;loadprops:multiloadproptype;
                           nlc,nb:integer;beamtop:beamtopoltype;Lv,Av,Iv,alfv:realbeamtype;E,alft:real;
                           xglob:reallcdoftype;npoints:integer;xs:realpointstype;
                           var defx,defy:reallcpointstype);

//to calculate x-y displacements in a given cross-section of a beam



var
   i,ip,ilc,j,k,n1,n2,nl:               integer;
   x,s,L,a,b,c,F,M,q1,q2,m1,m2:         real;
   x0,dxn,dyn,dxe,dye,dxl,dyl,dx,dy:    real;
   Ma,Mc,A0,C0:                         real;
   xloc,Txloc,fs:                       real6type;
   T,Kloc:                              stifloctype;

begin
   //preparation
   L:=Lv[nb];
   Transf(alfv[nb],T);

   //start dof, end dof
   n1:=(beamtop[nb,1]-1)*3;
   n2:=(beamtop[nb,2]-1)*3;

   for ilc:=1 to nlc do
      begin

         nl:=nloads[ilc,nb];

         //local displacement vector in global coordinate system
         for i:=1 to 3 do
            begin
               xloc[i]:=xglob[ilc,n1+i];
               xloc[i+3]:=xglob[ilc,n2+i];
            end;

         //transform displacements from global to local coordinate system
         for i:=1 to 6 do
            begin
               s:=0;
               for j:=1 to 6 do
                  s:=s+T[j,i]*xloc[j];    //T_transp ? xloc
               Txloc[i]:=s;
            end;
         for i:=1 to 6 do
            xloc[i]:=Txloc[i];

         //end moments from loads
         Ma:=0;
         Mc:=0;
         for j:=1 to nl do
            begin
               Calc_qloc(Lv[nb],Av[nb],E,alft,loadtyps[ilc,nb,j],loadprops[ilc,nb,j],fs);
               Ma:=Ma-fs[3];
               Mc:=Mc-fs[6];
            end;

         //cycle begins for x points
         for ip:=1 to npoints do
            begin
               x:=xs[ip];

               //x-y displacements from nodal displacments
               dxn:=xloc[1]*(1-x);
               dxn:=dxn+xloc[4]*x;
               dyn:=xloc[3]*Lv[nb]*( PowInt(x,3)-2*x*x+x );
               dyn:=dyn+xloc[6]*Lv[nb]*( PowInt(x,3)-x*x );
               dyn:=dyn+xloc[2]*( 2*PowInt(x,3)-3*x*x+1 );
               dyn:=dyn+xloc[5]*( -2*PowInt(x,3)+3*x*x );

               //y displacements from end moments (from loads on a clamped beam)
               dye:=PowInt(Lv[nb],2)/6*x*( Ma*(x*x-3*x+2) + Mc*(x*x-1) );
               dye:=dye/E/Iv[nb];
               dxe:=0;

               //x-y displacements from loads (but with zero end moments)
               x0:=x*L;
               dxl:=0;
               dyl:=0;
               for j:=1 to nl do
                  begin
                     case loadtyps[ilc,nb,j] of
                        1: begin  //conc force perp
                              F:=loadprops[ilc,nb,j,1];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=L-a;
                              if x0<=a then dyl:=dyl+F*c*x0*(a*a+2*a*c-x0*x0)/L/6;
                              if x0>a then dyl:=dyl-F*a*(PowInt(a,3)+a*a*c-3*x0*a*a+3*x0*x0*a-4*x0*a*c+3*x0*x0*c-PowInt(x0,3)-2*x0*c*c)/L/6;
                           end;
                        2: begin  //conc force parallel
                              F:=loadprops[ilc,nb,j,1];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=L-a;
                              if x0<=a then dxl:=dxl+c/L*F*x0;
                              if x0>a then dxl:=dxl+a/L*F*(L-x0);
                           end;
                        3: begin  //dist load perp
                              q1:=loadprops[ilc,nb,j,1];
                              q2:=loadprops[ilc,nb,j,2];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=(1-loadprops[ilc,nb,j,4])*L;
                              b:=L-a-c;
                              A0:=q2*b/2*(b/3+c)/L+q1*b/2*(2*b/3+c)/L;
                              C0:=q2*b/2*(a+2/3*b)/L+q1*b/2*(a+1/3*b)/L;
                              if x0<=a then
                                 dyl:=dyl+(1/2*A0*L*a*a-1/8*L*q1*b*b*b-1/24*L*q2*b*b*b+1/8*q1*a*b*b*b
                                        -2/3*C0*L*L*L+11/120*q1*b*b*b*b+1/30*q2*b*b*b*b+1/24*q2*a*b*b*b
                                        -1/6*A0*x0*x0*L-1/3*A0*a*a*a-1/3*A0*b*b*b-1/3*C0*a*a*a+L*L*C0*c
                                        -A0*b*a*a-A0*b*b*a+1/2*L*A0*b*b-C0*a*a*b-C0*a*b*b+C0*L*a*a
                                       +C0*L*b*b+L*A0*b*a+2*C0*L*a*b-1/3*C0*b*b*b)*x0/L;
                              if (x0>a) and (x0<(a+b)) then
                                 dyl:=dyl+1/120*(L*q2-L*q1)/b/L*PowInt(x0,5)
                                       +1/120*(5*L*q1*a+5*L*q1*b-5*L*q2*a)/b/L*PowInt(x0,4)
                                       +1/120*(-20*L*A0*b-10*L*q1*a*a-20*L*q1*a*b+10*L*q2*a*a)/b/L*PowInt(x0,3)
                                       +1/120*(10*L*q1*a*a*a-10*L*q2*a*a*a+30*L*q1*b*a*a)/b/L*x0*x0
                                       +1/120*(-5*L*q1*a*a*a*a-20*L*q1*a*a*a*b-5*L*q2*b*b*b*b+5*L*q2*a*a*a*a-120*A0*b*b*a*a
                                               -120*A0*b*b*b*a-40*A0*b*b*b*b-15*L*q1*b*b*b*b-40*A0*a*a*a*b+60*L*A0*b*a*a
                                               +120*L*A0*b*b*a+60*L*A0*b*b*b-40*b*b*b*b*C0+4*q2*b*b*b*b*b+11*q1*b*b*b*b*b
                                               -40*a*a*a*b*C0+5*q2*a*b*b*b*b-120*a*a*b*b*C0-120*a*b*b*b*C0+15*q1*a*b*b*b*b
                                               -80*C0*b*L*L*L+120*C0*b*b*b*L+120*L*L*C0*c*b+120*C0*b*L*a*a+240*C0*b*b*L*a)/b/L*x0
                                       +1/120*(L*q1*a*a*a*a*a-L*q2*a*a*a*a*a+5*L*q1*b*a*a*a*a)/b/L;
                              if x0>=(a+b) then
                                 dyl:=dyl+1/120*(L-x0)*(40*L*L*C0*x0-120*L*C0*a*b-60*L*C0*a*a
                                                     -60*L*C0*b*b-20*L*C0*x0*x0-5*q2*a*b*b*b
                                                     -4*q2*b*b*b*b+120*a*a*b*C0-11*q1*b*b*b*b
                                                     +120*A0*b*a*a+120*a*b*b*C0+40*C0*a*a*a
                                                     +40*b*b*b*C0-15*q1*a*b*b*b+40*A0*b*b*b
                                                     +40*a*a*a*A0+120*A0*b*b*a)/L;
                           end;
                        4: begin  //dist load parallel
                              q1:=loadprops[ilc,nb,j,1];
                              q2:=loadprops[ilc,nb,j,2];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=(1-loadprops[ilc,nb,j,4])*L;
                              b:=L-a-c;
                              A0:=1/6*b*(q2*b+2*q1*b+3*c*q1+3*c*q2)/L;
                              //C0:=s*(0.5*b+a)/L;
                              if x0<=a then
                                 dxl:=dxl-1/6*x0*(-6*A0*L+6*A0*a-q2*b*b-2*q1*b*b+6*A0*b+6*c*A0-3*b*q1*c-3*b*q2*c)/L;
                              if (x0>a) and (x0<(a+b)) then
                                 dxl:=dxl-1/6*(L*q2-L*q1)/b/L*x0*x0*x0
                                       -1/6*(-3*L*a*q2+3*L*q1*b+3*L*a*q1)/b/L*x0*x0
                                       -1/6*(3*a*a*q2*L-3*c*b*b*q2+6*c*b*A0-3*c*b*b*q1-q2*b*b*b
                                            +6*A0*a*b-3*a*a*q1*L+6*A0*b*b-2*q1*b*b*b-6*b*L*a*q1
                                            -6*A0*b*L)/b/L*x0
                                       -1/6*(-L*q2*a*a*a+L*q1*a*a*a+3*L*q1*b*a*a)/b/L;
                              if x0>=(a+b) then
                                 dxl:=dxl+1/6*(L-x0)*b*(3*a*q2+3*a*q1+2*b*q2+q1*b)/L;
                           end;
                        5: begin  //conc moment
                              M:=loadprops[ilc,nb,j,1];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=L-a;
                              if x0<=a then dyl:=dyl+M*x0*(x0*x0+3*a*a-4*L*L+6*c*L)/L/6;
                              if x0>a then dyl:=dyl-M*(L-x0)*(x0*x0-2*x0*L+3*a*a)/L/6;
                           end;
                        6: begin  //dist moment
                              m1:=loadprops[ilc,nb,j,1];
                              m2:=loadprops[ilc,nb,j,2];
                              a:=loadprops[ilc,nb,j,3]*L;
                              c:=(1-loadprops[ilc,nb,j,4])*L;
                              b:=L-a-c;
                              if x0<=a then
                                 dyl:=dyl+1/24*b*x0/L*(12*c*L*m1+12*c*L*m2+b*b*m1+3*b*b*m2+8*m1*b*L
                                                   +4*L*m2*b-8*L*L*m1-8*L*L*m2+2*m1*x0*x0
                                                    +2*m2*x0*x0+6*m1*a*a+6*m2*a*a+4*a*m1*b+8*m2*a*b);
                              if (x0>a) and (x0<(a+b)) then
                                 dyl:=dyl+1/24/L*(m1*L-m2*L)/b*x0*x0*x0*x0
                                       +1/24/L*(-4*m1*L*a+2*b*b*m2+4*m2*a*L-4*m1*b*L+2*b*b*m1)/b*x0*x0*x0
                                       +1/24/L*(6*a*a*L*m1-6*a*a*L*m2+12*a*L*m1*b)/b*x0*x0
                                       +1/24/L*(6*m1*a*a*b*b-4*m1*L*a*a*a+m1*b*b*b*b+6*m2*a*a*b*b
                                               +8*m2*a*b*b*b+4*m2*a*a*a*L-12*m1*b*L*a*a+12*b*b*c*L*m2
                                               -8*b*b*m1*L*L+3*m2*b*b*b*b+4*m1*a*b*b*b+8*m1*b*b*b*L
                                               -8*b*b*m2*L*L+4*L*m2*b*b*b+12*b*b*c*L*m1)/b*x0
                                       +1/24/L*(-m2*L*a*a*a*a+4*m1*b*L*a*a*a+m1*L*a*a*a*a)/b;
                              if x0>=(a+b) then
                                 dyl:=dyl-1/24*b*(L-x0)*(-4*x0*L*m2-4*x0*L*m1+b*b*m1+3*b*b*m2
                                               +6*m1*a*a+6*m2*a*a+4*a*m1*b+2*m2*x0*x0+8*m2*a*b+2*m1*x0*x0)/L;
                           end;
                        7: begin  //temperature change
                              //no displacements
                           end;
                     end; //case
                  end; //for j

               dxl:=dxl/E/Av[nb];
               dyl:=dyl/E/Iv[nb];

               dx:=dxn+dxe+dxl;
               dy:=dyn+dye+dyl;

               //tranform to global (???)
               dx:=dx;
               dy:=dy;

               //store dx,dy
               defx[ilc,ip]:=dx;
               defy[ilc,ip]:=dy;

            end; //for ip

      end; //for ilc

end; //CalcDeformations


function SetLMin(nbeam:integer;Lv:realbeamtype):real;

var
   i:      integer;
   Lmax:   real;

begin
   LMax:=0;
   for i:=1 to nbeam do
      Lmax:=MaxReal(Lmax,Lv[i]);
   SetLMin:=Lmax/LminDiv; //"18" helyett Options men?bol ?ll?tjuk
end; //SetLMin


procedure BeamEndOffset(nbeam:integer;beamtop:beamtopoltype;Lv,alfv:realbeamtype;
                        var offset:realbeam2type);

var
   ib,jb,iorig,jorig,nprof:   integer;
   ang_brace,ang_chord,h,bi:  real;
   chordconnect:   boolean;
   name:  string;

begin
   //zero out
   for ib:=1 to maxbeam do
      begin
         offset[ib,1]:=zero;
         offset[ib,2]:=1-zero;
      end;

   //set offsets
   for ib:=1 to nbeam do
      begin
         iorig:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
         if mainf.geom.ARud[iorig-1].tipus in [rBelsoOszlop,rKulsoOszlop,rRacs,rBak] then
            begin
               //inclination of the element (brace)
               ang_brace:=alfv[ib];

               //searching for adjoining chord members -- start node
               chordconnect:=False;
               jb:=0;
               repeat
                  jb:=jb+1;
                  if jb<>ib then
                     begin
                        if (beamtop[jb,1]=beamtop[ib,1]) or (beamtop[jb,2]=beamtop[ib,1]) then
                        //if (beamtop[jb,1]=beamtop[ib,2]) or (beamtop[jb,2]=beamtop[ib,2]) then
                           begin
                              jorig:=mainf.geom.aStatvaz[jb-1].RudIdx+1;
                              if (mainf.geom.ARud[jorig-1].tipus in [rFelsoOv,rAlsoOv]) and
                                 (not mainf.geom.aStatVaz[jb-1].Stagbar) then
                                 begin
                                    chordconnect:=True;

                                    //chord height
                                    name:=mainf.geom.aRud[jorig-1].profile.profile;
                                    nprof:=0;
                                    repeat
                                      nprof:=nprof+1;
                                    until aProfile[nprof-1].profile=name;
                                    h:=aProfile[nprof-1].H;

                                    //bolt distance
                                    if mainf.geom.aStatVaz[ib-1].PurlinSupport
                                       then bi:=mainf.geom.aRud[jorig-1].szelv[1].b1
                                       else bi:=mainf.geom.aRud[jorig-1].szelv[1].b2;

                                    //chord inclination
                                    ang_chord:=alfv[jb];

                                    //offset
                                    offset[ib,1]:=(h/2-bi)/Abs(Sin(ang_chord-ang_brace))/Lv[ib];
                                 end;
                           end;
                     end;
               until (jb=nbeam) or chordconnect;

               //searching for adjoining chord members -- end node
               chordconnect:=False;
               jb:=0;
               repeat
                  jb:=jb+1;
                  if jb<>ib then
                     begin
                        //if (beamtop[jb,1]=beamtop[ib,1]) or (beamtop[jb,2]=beamtop[ib,1]) then
                        if (beamtop[jb,1]=beamtop[ib,2]) or (beamtop[jb,2]=beamtop[ib,2]) then
                           begin
                              jorig:=mainf.geom.aStatvaz[jb-1].RudIdx+1;
                              if (mainf.geom.ARud[jorig-1].tipus in [rFelsoOv,rAlsoOv]) and
                                 (not mainf.geom.aStatVaz[jb-1].Stagbar) then
                                 begin
                                    chordconnect:=True;

                                    //chord height
                                    name:=mainf.geom.aRud[jorig-1].profile.profile;
                                    nprof:=0;
                                    repeat
                                      nprof:=nprof+1;
                                    until aProfile[nprof-1].profile=name;
                                    h:=aProfile[nprof-1].H;

                                    //bolt distance
                                    if mainf.geom.aStatVaz[ib-1].PurlinSupport
                                       then bi:=mainf.geom.aRud[jorig-1].szelv[1].b1
                                       else bi:=mainf.geom.aRud[jorig-1].szelv[1].b2;

                                    //chord inclination
                                    ang_chord:=alfv[jb];

                                    //offset
                                    offset[ib,2]:=1-(h/2-bi)/Abs(Sin(ang_chord-ang_brace))/Lv[ib];
                                 end;
                           end;
                     end;
               until (jb=nbeam) or chordconnect;

            end;

      end; //for ib

end; //BeamEndOffset


procedure BeamCheckPoints(nloads:intlcbeamtype;loadtyps:intlcbeamloadtype;
                          loadprops:multiloadproptype;nlc,nb:integer;L,Lmin:real;
                          offset:real2type;var np:integer;var xrel:realpointstype);

//to define the position of checkpoints within a beam

var
   ip,ilc,i,j,nl,nch,nadd:  integer;
   dL:                      real;

begin
   Lmin:=Lmin/L;

   ip:=1;
   xrel[ip]:=zero;
   ip:=ip+1;
   xrel[ip]:=offset[1];
   ip:=ip+1;
   xrel[ip]:=0.5;

   //considering loads
   for ilc:=1 to nlc do
      begin
         nl:=nloads[ilc,nb];
         for j:=1 to nl do
            begin
               case loadtyps[ilc,nb,j] of
                  1,2,5: begin  //conc force or moment
                            if (loadprops[ilc,nb,j,3]-zero)>0 then
                               begin
                                  ip:=ip+1;
                                  xrel[ip]:=loadprops[ilc,nb,j,3]-zero;
                               end;
                            if (loadprops[ilc,nb,j,3]+zero)<1 then
                               begin
                                  ip:=ip+1;
                                  xrel[ip]:=loadprops[ilc,nb,j,3]+zero;
                               end;
                         end;
                  3,4,6: begin  //dist load or moment
                            if ((loadprops[ilc,nb,j,3]-zero)>0) and ((loadprops[ilc,nb,j,3]+zero)<1) then
                               begin
                                  ip:=ip+1;
                                  xrel[ip]:=loadprops[ilc,nb,j,3];
                               end;
                            if ((loadprops[ilc,nb,j,4]-zero)>0) and ((loadprops[ilc,nb,j,4]+zero)<1) then
                               begin
                                  ip:=ip+1;
                                  xrel[ip]:=loadprops[ilc,nb,j,4];
                               end;
                         end;
                 7: begin  //temperature change
                        //no checkpoints are added
                    end;
               end; //case
            end; //for j
      end; //for ilc

   ip:=ip+1;
   xrel[ip]:=offset[2];
   ip:=ip+1;
   xrel[ip]:=1-zero;

   //ordering
   repeat
      nch:=0;
      for j:=1 to (ip-1) do
         if xrel[j]>xrel[j+1] then
            begin
               ChangeReal(xrel[j],xrel[j+1]);
               nch:=nch+1;
            end;
   until nch=0;

   //add more points
   j:=1;
   while j<ip do
      begin
         dL:=xrel[j+1]-xrel[j];
         if dL>Lmin then
            begin
               nadd:=Trunc(dL/Lmin);
               for i:=ip downto (j+1) do
                  xrel[i+nadd]:=xrel[i];
               for i:=1 to nadd do
                  xrel[j+i]:=xrel[j]+i*dL/(nadd+1);
               ip:=ip+nadd;
               j:=j+nadd;
            end; //if
         j:=j+1;
      end;  // while

   np:=ip;

end; //BeamCheckPoints



procedure ReduceCheckPoints(nlc:integer;L:real;offset:real2type;var np:integer;
                            var xs:realpointstype;var Nel,Tel,Mel:reallcpointstype;
                            var activecp:int2type);

//this routine eliminates the duplicated checkpoints

var
   ip,il,j:   integer;
   elim:  boolean;

begin
   ip:=2;
   while ip<np do  //possible elimination of points 2 to (np-1)
      begin
         elim:=True;
         il:=0;
         repeat
            il:=il+1;
            if
               ( (Abs(Nel[il,ip]-Nel[il,ip-1])>0.001) or    //compare to previous point
                 (Abs(Tel[il,ip]-Tel[il,ip-1])>0.001) or
                 (Abs(Mel[il,ip]-Mel[il,ip-1])>0.001) )
               and
               ( (Abs(Nel[il,ip]-Nel[il,ip+1])>0.001) or    //compare to next point
                 (Abs(Tel[il,ip]-Tel[il,ip+1])>0.001) or
                 (Abs(Mel[il,ip]-Mel[il,ip+1])>0.001) )
            then elim:=False;
         until (il=nlc) or (not elim);  //il for load cases
         if elim then
            begin
               np:=np-1;
               for j:=ip to np do
                  begin
                     xs[j]:=xs[j+1];
                     for il:=1 to nlc do
                        begin
                           Nel[il,j]:=Nel[il,j+1];
                           Tel[il,j]:=Tel[il,j+1];
                           Mel[il,j]:=Mel[il,j+1];
                        end;
                  end;
               ip:=ip-1;
            end;
         ip:=ip+1;
      end; //ip for points

   activecp[1]:=1;
   activecp[2]:=np;
   for ip:=1 to np do
      if xs[ip]<=offset[1] then activecp[1]:=ip;
   for ip:=np downto 1 do
      if xs[ip]>=offset[2] then activecp[2]:=ip;

end; //ReduceCheckPoints




procedure StatDocShort(nlc,nbeam,nnode:integer;npoints:intbeamtype;xrel:realbeampointstype;
                       beamtop:beamtopoltype;alfv:realbeamtype;name:shortstring;
                       Nel,Tel,Mel,defx,defy:realbeamlcpointstype);

var
   i,j,j1,j2,ilc,ib,icp,irow,inode,ct,cmax:  integer;
   Nmin,Nmax,Tmin,Tmax,Mmin,Mmax,dxmin,dxmax,dymin,dymax,dx,dy:  real;
   umax,xmax,ULSmax,SLSmax,Rxmax,Rxmin,Rymax,Rymin:  real;
   f: text;
   s,smax: string;
   v: string10type;
   ib2: array [1..8] of integer;
   side: Toldal;


begin
   Assign(f,name);
   ReWrite(f);


   //title
   i:=0;
   j:=0;
   repeat
      i:=i+1;
      if name[i]='\' then j:=i;
   until i=Length(name);
   s:=Copy(name,j+1,i-j);
   j1:=Pos('_brief',s);
   s:=Copy(s,1,j1-1);
   WriteLn(f,'**************************************************************');
   WriteLn(f,' ************************************************************');
   WriteLn(f);
   WriteLn(f,' BRIEF DOCUMENTATION FOR THE STRUCTURAL DESIGN OF LINDABTRUSS');
   WriteLn(f);
   for i:=1 to ((62-Length(s)) div 2) do
      Write(f,' ');
   WriteLn(f,s);
   WriteLn(f,'                         ',DateTimeToStr(date));
   WriteLn(f);
   WriteLn(f,' ************************************************************');
   WriteLn(f,'**************************************************************');
   WriteLn(f);
   WriteLn(f);


   //input data

   WriteLn(f,'**********');
   WriteLn(f);
   WriteLn(f,'INPUT DATA');
   WriteLn(f);
   WriteLn(f,'**********');
   WriteLn(f);

   with mainf do
      begin
         WriteLn(f,'*** Main dimensions');
         WriteLn(f,'  Span:                          ',FloatToStr(EditL.value),' mm');
         if not r_calceave.checked then WriteLn(f,'  Eave height:                   ',FloatToStr(EditVall.value),' mm')
                           else WriteLn(f,'  Eave height:                   ','calculated');
         if not r_calcridge.checked then WriteLn(f,'  Ridge height:                  ',FloatToStr(EditTarej.value),' mm')
                                    else WriteLn(f,'  Ridge height:                  ','calculated');
         if not r_calcSlope.checked
            then begin
                    if BtnSlopePerc.Down then WriteLn(f,'  Slope:                         ',FloatToStr(Round(TanD(EditSlopeDeg.value)*10000)/100),' %');
                    if BtnSlopeDeg.Down then WriteLn(f,'  Slope:                         ',FloatToStr(EditSlopeDeg.value),' deg');
                 end
            else WriteLn(f,'  Slope:                         ','calculated');
         if c_no_eaves.checked then WriteLn(f,'  Eave type:                     ','no eaves');
         if c_stag_bar.checked then WriteLn(f,'  Eave type and size:            ','eave console:  ',FloatToStr(EditEave.value),' mm');
         if c_overhang.checked then WriteLn(f,'  Eave type and size:            ','overhang:  ',FloatToStr(EditEave.value),' mm');
         if c_gorgo_jobb.Checked then WriteLn(f,'  Supports:                      ','hinge on the left, roller on the right');
         if c_gorgo_bal.Checked then WriteLn(f,'  Supports:                       ','roller on the left, hinge on the right');
         if c_nincs_gorgo.Checked then WriteLn(f,'  Supports:                     ','hinges at both supports');
         WriteLn(f,'  Purlin distribution:           ',FloatToStr(EditDist1.value),' mm + ',FloatToStr(EditN.value),' ? ',FloatToStr(EditDist.value),' mm + ',FloatToStr(EditDist2.value),' mm');
         if c_Horizontal.Checked then WriteLn(f,'  Purlin positions are measured: ','horizontally');
         if c_on_roof.Checked then WriteLn(f,'  Purlin positions are measured: ','on the roof');
         WriteLn(f);

         WriteLn(f,'*** Chords, brace members, bolts');
         WriteLn(f,'  Basic upper chord profile:                           ',c_def_fov.text);
         WriteLn(f,'  Upper chord bolt distances (from outside to inside): ',FloatToStr(SpinTopB1.value),' mm + ',FloatToStr(SpinTopRows.value-1),' ? ',FloatToStr(EditTopDist.value),' mm + ',FloatToStr(SpinTopB2.value),' mm');
         if c_upper_purlins.checked then
            begin
               WriteLn(f,'  Lateral supports along upper chord:                  ','at purlins');
               if chkBak.checked then WriteLn(f,'  Purlin support size:                                 ',FloatToStr(EditBak.value),' mm');
            end;
         if c_upper_cont.checked then
            begin
               WriteLn(f,'  Lateral supports along upper chord:                  ','continuous');
               if chkBak.checked then WriteLn(f,'  Load distance at upper chord:                        ',FloatToStr(EditDistUpper.value),' mm');
            end;

         WriteLn(f,'  Basic lower chord profile:                           ',c_def_aov.text);
         WriteLn(f,'  Lower chord bolt distances (from outside to inside): ',FloatToStr(SpinBottomB1.value),' mm + ',FloatToStr(SpinBottomRows.value-1),' ? ',FloatToStr(EditBottomDist.value),' mm + ',FloatToStr(SpinBottomB2.value),' mm');
         if c_lower_cont.checked then WriteLn(f,'  Lateral supports along lower chord:                  ','continuous');
         if c_lower_each.checked then WriteLn(f,'  Lateral supports along lower chord:                  ','each node');
         if c_lower_none.checked then WriteLn(f,'  Lateral supports along lower chord:                  ','none');
         if c_lower_custom.checked then WriteLn(f,'  Lateral supports along lower chord:                ','user defined');
         WriteLn(f,'  Load distance at lower chord:                        ',FloatToStr(EditDistLower.value),' mm');

         WriteLn(f,'  Basic brace member profile:                          ',c_def_racs.text);
         if chkMer.checked then Write(f,'  Columns:                                             perpendicular to upper chord, ');
         if chkFugg.checked then Write(f,'  Columns:                                             vertical, ');
         if StrToInt(EditEvery.text)=1 then Writeln(f,'at every purlin');
         if StrToInt(EditEvery.text)=2 then Writeln(f,'at every 2nd purlin');
         if StrToInt(EditEvery.text)=3 then Writeln(f,'at every3rd purlin');
         if StrToInt(EditEvery.text)=4 then Writeln(f,'at every 4th purlin');
         Write(f,'  Column shift:                                        ');
         if StrToInt(EditSkip.text)=0 then WriteLn(f,'no shift');
         if StrToInt(EditSkip.text)=1 then WriteLn(f,'shift by ',EditSkip.text,' purlin');
         if StrToInt(EditSkip.text)>1 then WriteLn(f,'shift by ',EditSkip.text,' purlins');

         WriteLn(f,'  Bolt diamater and grade:                             ',c_bolts_dia.text,', ',c_bolts_grade.text);

         WriteLn(f);
         WriteLn(f,'*** Profiles');
         WriteLn(f,'  See the drawings.');

         WriteLn(f);
         WriteLn(f,'*** Load cases');
         for ilc:=1 to nlc do
            begin
       { TODO : fagy 4-n?l }        WriteLn(f,'  Load case: ',geom.aLoad[ilc-1].LC,' in load group: ',c_groupmod.items[ c_GroupMod.items.indexof(  geom.LCGroups[ilc-1] )] );

               for irow:=0 to geom.aLoad[ilc-1].lines.count-1 do
                  begin
                     s:=geom.aLoad[ilc-1].lines[irow];
                     StringIn3(s,v);
                     case StrToInt(v[1]) of
                        1: WriteLn(f,'    type: 1, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m');
                        2: WriteLn(f,'    type: 2, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m');
                        3: WriteLn(f,'    type: 3, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m');
                        4: WriteLn(f,'    type: 4, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m',', angle: ',StrToFloat(v[9]):5:2,' deg');
                        8: WriteLn(f,'    type: 8, ','acting at: ',v[2]:20,',  F: ',StrToFloat(v[3]):5:1,' kN/m',', angle: ',StrToFloat(v[9]):5:2,' deg',', pos: ',v[10]);
                     end; //case
                  end;
            end;

         WriteLn(f);
         WriteLn(f,'*** Load combinations');
         WriteLn(f,'  Combinations are in accordance with: ',c_env.text);
         WriteLn(f,'  The applied partial and combination factors as follows:');
         for j:=1 to groups.ColCount do
            begin
               Write(f,'  ',STrunc(Groups.cells[0,j],30):30,STrunc(Groups.cells[1,j],15):15);
               if groups.Cells[1,j]=get('Permanent') then
                  begin
                     Write(f,'  ',' gam_G,inf=',Groups.cells[2,j]);
                     Write(f,'  ',' zeta =',Groups.cells[3,j]);
                     Write(f,'  ',' gam_G,sup=',Groups.cells[4,j]);
                     Write(f,'  ',' gam_G,A=',Groups.cells[5,j]);
                //     for i:=2 to 5 do
                //        Write(f,'  ',Groups.cells[i,j]);
                     WriteLn(f);
                  end;
               if groups.Cells[1,j]=get('Variable') then
                  begin
                     Write(f,'  ',' gam_Q    =',Groups.cells[7,j]);
                     Write(f,'  ',' psi_0=',Groups.cells[8,j]);
                     Write(f,'  ',' psi_1    =',Groups.cells[9,j]);
                     Write(f,'  ',' psi_2  =',Groups.cells[10,j]);
                //     for i:=7 to 10 do
                //        Write(f,'  ',Groups.cells[i,j]);
                     WriteLn(f);
                  end;
               if groups.Cells[1,j]=get('Accidental') then
                  begin
                   //  WriteLn(f,'  ',Groups.cells[6,j]);
                     WriteLn(f,'  ',' gam_A=    ',Groups.cells[6,j]);
                  end;
            end;
         WriteLn(f);


         WriteLn(f);
         WriteLn(f,'*** Checkpoints');
         //WriteLn(f,'  Checkpoint density used in calculation: ',AnsiReplaceStr( TMenuitem( Average1 ).Caption, '&',''));
         Write(f,'  Checkpoint density used in calculation: ');
         if StrBegin(BtnCheckPoints.Caption,Get('Checkpoints: Rare')) then WriteLn(f,'rare');
         if BtnCheckPoints.Caption=Get('Checkpoints: Average') then WriteLn(f,'average');
         if BtnCheckPoints.Caption=Get('Checkpoints: Dense') then WriteLn(f,'dense');

      end; //with

   WriteLn(f);
   WriteLn(f);
   WriteLn(f);


   //load combinations

   WriteLn(f,'************************');
   WriteLn(f);
   WriteLn(f,'ACTIVE LOAD COMBINATIONS');
   WriteLn(f);
   WriteLn(f,'************************');
   WriteLn(f);

   WriteLn(f);
   WriteLn(f,'*** ULS combinations');
   WriteLn(f);
   //WriteLn(f,'name   ','      permanent load',' + variable loads ');

   for ilc:=1 to mainf.comb.rowcount-1 do
      if (mainf.comb.cells[mainf.comb.colcount-1,ilc]='On') and
         (mainf.comb.cells[mainf.comb.colcount-2,ilc]='ULS') then
         begin
            Write(f,'  ',STrunc(mainf.comb.cells[0,ilc],10):10);
            for i:=1 to mainf.comb.colcount-3 do
               begin
                  if i=1 then Write(f,'   ',mainf.comb.cells[i,ilc])
                         else begin
                                 if (i mod 2)=1 then Write(f,' + ',mainf.comb.cells[i,ilc])
                                                else Write(f,'?','"',STrunc(mainf.comb.cells[i,ilc],10):10,'"');
                              end;
               end;
            WriteLn(f);
         end;

   WriteLn(f);
   WriteLn(f,'*** SLS combinations');
   WriteLn(f);
   //WriteLn(f,'name   ','      permanent load',' + variable loads ');

   for ilc:=1 to mainf.comb.rowcount-1 do
      if (mainf.comb.cells[mainf.comb.colcount-1,ilc]='On') and
         (mainf.comb.cells[mainf.comb.colcount-2,ilc]='SLS') then
         begin
            Write(f,'  ',STrunc(mainf.comb.cells[0,ilc],10):10);
            for i:=1 to mainf.comb.colcount-3 do
               begin
                  if i=1 then Write(f,'   ',mainf.comb.cells[i,ilc])
                         else begin
                                 if (i mod 2)=1 then Write(f,' + ',mainf.comb.cells[i,ilc])
                                                else Write(f,'?','"',STrunc(mainf.comb.cells[i,ilc],10):10,'"');
                              end;
               end;
            WriteLn(f);
         end;

   WriteLn(f);
   WriteLn(f);


   //internal forces, displacements

   WriteLn(f,'******************************************************');
   WriteLn(f);
   WriteLn(f,'MAX/MIN INTERNAL FORCES, DISPLACEMENTS FROM LOAD CASES');
   WriteLn(f);
   WriteLn(f,'******************************************************');
   WriteLn(f);
   WriteLn(f,'elem','    Nmin','    Nmax','    Vmin','    Vmax','    Mmin','    Mmax','   dxmin','   dxmax','   dymin','   dymax');
   WriteLn(f,'   #','    [kN]','    [kN]','    [kN]','    [kN]','   [kNm]','   [kNm]','    [mm]','    [mm]','    [mm]','    [mm]');
   WriteLn(f);

   for ib:=1 to nbeam do
      begin
         Nmin:=Nel[ib,1,1];
         Nmax:=Nel[ib,1,1];
         Tmin:=Tel[ib,1,1];
         Tmax:=Tel[ib,1,1];
         Mmin:=Mel[ib,1,1];
         Mmax:=Mel[ib,1,1];
         dx:=defx[ib,1,1];
         dy:=defy[ib,1,1];
         LtoG(dx,dy,alfv[ib]);
         dxmin:=dx;
         dxmax:=dx;
         dymin:=dy;
         dymax:=dy;
         for ilc:=1 to nlc do
            begin
               for icp:=1 to npoints[ib] do
                  begin
                     if Nel[ib,ilc,icp]<Nmin then Nmin:=Nel[ib,ilc,icp];
                     if Nel[ib,ilc,icp]>Nmax then Nmax:=Nel[ib,ilc,icp];
                     if Tel[ib,ilc,icp]<Tmin then Tmin:=Tel[ib,ilc,icp];
                     if Tel[ib,ilc,icp]>Tmax then Tmax:=Tel[ib,ilc,icp];
                     if Mel[ib,ilc,icp]<Mmin then Mmin:=Mel[ib,ilc,icp];
                     if Mel[ib,ilc,icp]>Mmax then Mmax:=Mel[ib,ilc,icp];
                     dx:=defx[ib,ilc,icp];
                     dy:=defy[ib,ilc,icp];
                     LtoG(dx,dy,alfv[ib]);
                     if dx<dxmin then dxmin:=dx;
                     if dx>dxmax then dxmax:=dx;
                     if dy<dymin then dymin:=dy;
                     if dy>dymax then dymax:=dy;
                  end; //checkpoint
            end; //load case

            WriteLn(f,ib:4,Nmin/1000:8:2,Nmax/1000:8:2,Tmin/1000:8:2,Tmax/1000:8:2,Mmin/1000000:8:2,Mmax/1000000:8:2,dxmin:8:2,dxmax:8:2,dymin:8:2,dymax:8:2);

      end;  //beam

   WriteLn(f);
   WriteLn(f);

   
   //verification results

   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'MAXIMAL UTILIZATIONS IN ULS');
   WriteLn(f);
   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'elem #','   utilization','     combination','   position','    governing check');
   WriteLn(f);

   ULSmax:=-1;
   for ib:=1 to nbeam do
      begin
         umax:=ell[ib,1,1,1].util;
         smax:=ell[ib,1,1,1].tonkr;
         cmax:=KombKezd;
         xmax:=xrel[ib,1];
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  for icp:=1 to npoints[ib] do
                     begin
                        for i:=1 to length(ell[ib,ilc,icp]) do
                        if ell[ib,ilc,icp,i].util>umax then
                           begin
                              umax:=ell[ib,ilc,icp,i].util;
                              smax:=ell[ib,ilc,icp,i].tonkr;
                              cmax:=ilc;
                              xmax:=xrel[ib,icp];
                           end;

                     end; //checkpoint
               end; //load combination
               if umax>ULSmax then ULSmax:=umax;

               if Abs(umax)<0.00001
                  then WriteLn(f,ib:5,'         --    ','    --             ','    --   ','   not checked')
                  else WriteLn(f,ib:5,' ',umax*100:10:2,' %      ',STrunc(aLCorComb[cmax].name,10):10,' ',xmax*100:10:2,' %    ',smax)

         end;  //beam

   WriteLn(f);
   WriteLn(f,'Maximal ULS utilization:  ',ULSmax*100:10:2,' %');

   WriteLn(f);
   WriteLn(f);

   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'MAXIMAL UTILIZATIONS IN SLS');
   WriteLn(f);
   WriteLn(f,'***************************');
   WriteLn(f);

   WriteLn(f,'Displacement limit:  L/',mainf.DeflLimit.text);
   WriteLn(f);

   SLSMax:=0;
   for ilc:=KombKezd to KombVeg do
      if aLCOrComb[ilc].sls then
         begin
            if aLCorComb[ilc].MaxSLSChkRatio>SLSMax then SLSMax:=aLCorComb[ilc].MaxSLSChkRatio;
            WriteLn(f,STrunc(aLCorComb[ilc].name,10):10,': ',aLCorComb[ilc].MaxSLSChkRatio:10:2,' %');
         end;
   WriteLn(f);
   WriteLn(f,'Maximal SLS utilization:  ',SLSmax:10:2,' %');
   WriteLn(f);
   WriteLn(f);


   //output for extr checks

   WriteLn(f,'******************************');
   WriteLn(f);
   WriteLn(f,'OUTPUT DATA FOR FURTHER CHECKS');
   WriteLn(f);
   WriteLn(f,'******************************');
   WriteLn(f);

   with mainf do
      begin
         WriteLn(f,'*** Max/min reaction forces');
         for ib:=1 to 2 do
            begin
               Rxmax:=0;
               Rxmin:=0;
               Rymax:=0;
               Rymin:=0;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if reac[ilc,ib,1]>Rxmax then Rxmax:=reac[ilc,ib,1];
                        if reac[ilc,ib,1]<Rxmin then Rxmin:=reac[ilc,ib,1];
                        if reac[ilc,ib,2]>Rymax then Rymax:=reac[ilc,ib,2];
                        if reac[ilc,ib,2]<Rymin then Rymin:=reac[ilc,ib,2];
                     end;
               if ib=1 then Write(f,'  Left support,  ')
                       else Write(f,'  Right support, ');
               WriteLn(f,'horizontal reactions:   ','Rxmin=',Rxmin/1000:8:1,' kN    ','Rxmax=',Rxmax/1000:8:1,' kN');
               if ib=1 then Write(f,'  Left support,  ')
                       else Write(f,'  Right support, ');
               WriteLn(f,'vertical reactions:     ','Rymin=',Rymin/1000:8:1,' kN    ','Rymax=',Rymax/1000:8:1,' kN');
            end;
         WriteLn(f);


         WriteLn(f,'*** Internal forces on lower chord connection');

         //finding lower chord chord connection
         icp:=-1;
         inode:=0;
         repeat
            inode:=inode+1;
            ct:=0;
            for ib:=1 to nbeam do
               if (beamtop[ib,1]=inode) or (beamtop[ib,2]=inode) then
                  begin
                     ct:=ct+1;
                     ib2[ct]:=ib;
                  end;
            if (ct=2) and (not geom.aStatCsp[inode-1].CsuklosTam) and (not geom.aStatCsp[inode-1].GorgosTam) then
               begin
                  ib:=ib2[1];
                  j1:=mainf.geom.aStatvaz[ib2[1]-1].RudIdx+1;
                  j2:=mainf.geom.aStatvaz[ib2[2]-1].RudIdx+1;
                  if (mainf.geom.ARud[j1-1].tipus in [rAlsoOv]) and (mainf.geom.ARud[j2-1].tipus in [rAlsoOv]) then
                     begin
                        if beamtop[ib,1]=inode then icp:=1;
                        if beamtop[ib,2]=inode then icp:=npoints[ib];
                      end;
               end;
         until (icp>0) or (i=nnode);

         //internal forces: node=inode, beam=ib, checkpoint=icp

         //max N
         Nmax:=-1000000000;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Nel[ib,ilc,icp]>Nmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Max N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
         //min N
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Nel[ib,ilc,icp]<Nmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Min N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

         //max T
         Tmax:=-1000000000;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Tel[ib,ilc,icp]>Tmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Max V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
         //min T
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Tel[ib,ilc,icp]<Tmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Min V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

         //max M
         Mmax:=-1000000000;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Mel[ib,ilc,icp]>Mmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Max M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
         //min M
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Mel[ib,ilc,icp]<Mmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Min M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

         WriteLn(f);


         WriteLn(f,'*** Internal forces on ridge connection');

         //finding ridge
         icp:=-1;
         inode:=0;
         repeat
            inode:=inode+1;
            ct:=0;
            for ib:=1 to nbeam do
               if (beamtop[ib,1]=inode) or (beamtop[ib,2]=inode) then
                  begin
                     ct:=ct+1;
                     ib2[ct]:=ib;
                  end;
            if (ct=2) and (mainf.geom.aStatVaz[ib2[1]-1].oldal<>mainf.geom.aStatVaz[ib2[2]-1].oldal) then
               begin
                  j1:=mainf.geom.aStatvaz[ib2[1]-1].RudIdx+1;
                  j2:=mainf.geom.aStatvaz[ib2[2]-1].RudIdx+1;
                  if (mainf.geom.ARud[j1-1].tipus in [rFelsoOv]) and (mainf.geom.ARud[j2-1].tipus in [rFelsoOv]) then icp:=1;
               end;
         until (icp>0) or (i=nnode);


         for i:=1 to 2 do
            begin
               ib:=ib2[i];
               if beamtop[ib,1]=inode then icp:=1;
               if beamtop[ib,2]=inode then icp:=npoints[ib];

               //internal forces: node=inode, beam=ib, checkpoint=icp

               if mainf.geom.aStatVaz[ib-1].oldal=oBal then WriteLn(f,'  Internal forces from the left');
               if mainf.geom.aStatVaz[ib-1].oldal=oJobb then WriteLn(f,'  Internal forces from the right');

               //max N
               Nmax:=-1000000000;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Nel[ib,ilc,icp]>Nmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Max N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
               //min N
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Nel[ib,ilc,icp]<Nmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Min N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

               //max T
               Tmax:=-1000000000;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Tel[ib,ilc,icp]>Tmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Max V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
               //min T
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Tel[ib,ilc,icp]<Tmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Min V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

               //max M
               Mmax:=-1000000000;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Mel[ib,ilc,icp]>Mmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Max M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
               //min M
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Mel[ib,ilc,icp]<Mmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Min M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

               WriteLn(f);
            end; //for i:=1 to 2


         WriteLn(f,'*** Internal forces on eaves');

         if c_no_eaves.checked
            then
               begin
                  WriteLn(f,'  There are no eaves.');
               end

         else
            begin
               for i:=1 to 2 do
                  begin
                     case i of
                        1: side:=oBal;
                        2: side:=oJobb;
                     end;

                     //finding left/right eave
                     if c_overhang.checked then     //overhang
                        begin
                           ib:=0;
                           repeat
                              ib:=ib+1;
                              j1:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                           until (geom.ARud[j1-1].tipus=rFelsoov) and (geom.aStatVaz[ib-1].oldal=side)
                              and (geom.aStatVaz[ib-1].Overhang) and (geom.aStatVaz[ib-1].tipus=rFelsoOv)
                              and (CheckParallel(geom.ARud[j1-1].l,geom.aStatVaz[ib-1].v)) or (ib=nbeam);
                        end; //overhang
                     if c_stag_bar.checked then     //stag bar
                        begin
                           ib:=0;
                           repeat
                              ib:=ib+1;
                              j1:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                           until (geom.ARud[j1-1].tipus=rFelsoov) and (geom.aStatVaz[ib-1].oldal=side)
                              and (geom.aStatVaz[ib-1].Stagbar) and (geom.aStatVaz[ib-1].tipus=rFelsoOv)
                              and (CheckParallel(geom.ARud[j1-1].l,geom.aStatVaz[ib-1].v)) or (ib=nbeam);
                        end; //overhang)

                     case i of
                        1: j2:=1;
                        2: j2:=-1;
                     end;
                     if (j2*mainf.geom.aStatVaz[ib-1].v.x1>j2*mainf.geom.aStatVaz[ib-1].v.x2)
                        then icp:=1
                        else icp:=npoints[ib];


                     //internal forces: beam=ib, checkpoint=icp

                     case i of
                        1: WriteLn(f,'  Internal forces on the the left eave');
                        2: WriteLn(f,'  Internal forces on the the right eave');
                     end;

                     //max N
                     Nmax:=-1000000000;
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Nel[ib,ilc,icp]>Nmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Max N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
                     //min N
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Nel[ib,ilc,icp]<Nmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Min N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

                     //max T
                     Tmax:=-1000000000;
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Tel[ib,ilc,icp]>Tmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Max V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
                     //min T
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Tel[ib,ilc,icp]<Tmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Min V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

                     //max M
                     Mmax:=-1000000000;
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Mel[ib,ilc,icp]>Mmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Max M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
                     //min M
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Mel[ib,ilc,icp]<Mmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Min M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

                     WriteLn(f);
                  end; //for i:=1 to 2

            end; //if there are eaves

      end;  //with


   Close(f);
end;  //StatDocShort



procedure StatDocFull(nlc,nbeam,nnode:integer;npoints:intbeamtype;xrel:realbeampointstype;
                      beamtop:beamtopoltype;alfv:realbeamtype;name:shortstring;
                      Nel,Tel,Mel,defx,defy:realbeamlcpointstype);

var
   i,j,j1,j2,ilc,ib,jb,icp,irow,inode,ct,cmax:  integer;
   Nmin,Nmax,Tmin,Tmax,Mmin,Mmax,dxmin,dxmax,dymin,dymax,dx,dy:  real;
   umax,xmax,ULSmax,SLSmax,Rxmax,Rxmin,Rymax,Rymin:  real;
   f: text;
   s,smax: string;
   v: string10type;
   ib2: array [1..8] of integer;
   nonstruct,braceconnect,yet,yet2: boolean;
   side: Toldal;


function CheckConv(s:shortstring):integer;

   var
      c: integer;

   begin
      c:=0;
      if StrBegin(s,'brace, pure axial compression') then c:=1;
      if StrBegin(s,'brace, axial compression and bending - 1') then c:=2;
      if StrBegin(s,'brace, axial compression and bending - 2') then c:=3;
      if StrBegin(s,'brace, axial compression and bending - 3') then c:=4;
      if StrBegin(s,'brace, axial compression and bending - 4') then c:=5;
      if StrBegin(s,'brace, flexural buckling') then c:=6;
      if StrBegin(s,'brace, torsional or torsional-flexural buckling') then c:=7;
      if StrBegin(s,'brace, interaction of flexural buckling and bending - 1') then c:=8;
      if StrBegin(s,'brace, interaction of flexural buckling and bending - 2') then c:=9;
      if StrBegin(s,'brace, tension - plastic resistance') then c:=10;
      if StrBegin(s,'brace, tension - ultimate resistance - 1') then c:=11;
      if StrBegin(s,'brace, tension - ultimate resistance - 2') then c:=12;
      if StrBegin(s,'brace, axial tension and bending - 1') then c:=13;
      if StrBegin(s,'brace, axial tension and bending - 2') then c:=14;
      if StrBegin(s,'brace, axial tension and bending - 3') then c:=15;
      if StrBegin(s,'brace, axial tension and bending - 4') then c:=16;
      if StrBegin(s,'chord, pure axial compression') then c:=17;
      if StrBegin(s,'chord, interaction of shear and axial compression') then c:=18;
      if StrBegin(s,'chord, axial compression and bending - 1') then c:=19;
      if StrBegin(s,'chord, axial compression and bending - 2') then c:=20;
      if StrBegin(s,'chord, flexural buckling') then c:=21;
      if StrBegin(s,'chord, torsional or torsional-flexural buckling') then c:=22;
      if StrBegin(s,'chord, lateral-torsinal buckling') then c:=23;
      if StrBegin(s,'chord, interaction of flexural buckling and lateral-torsional buckling') then c:=24;
      if StrBegin(s,'lower chord, global buckling of built-up CS') then c:=25;
      if StrBegin(s,'chord, plastic resistance for pure axial tension') then c:=26;
      if StrBegin(s,'chord, interaction of shear and axial tension') then c:=27;
      if StrBegin(s,'chord, axial tension and bending - 1') then c:=28;
      if StrBegin(s,'chord, axial tension and bending - 2') then c:=29;
      if StrBegin(s,'upper chord, N+M, restrained flange') then c:=30;
      if StrBegin(s,'upper chord, N+M, free flange') then c:=31;
      if StrBegin(s,'upper chord, buckling of the free flange') then c:=32;
      if StrBegin(s,'upper chord, interaction of in-plane flexural buckling and bending') then c:=33;
      if StrBegin(s,'lower chord, N+M, restrained flange') then c:=34;
      if StrBegin(s,'lower chord, N+M, free flange') then c:=35;
      if StrBegin(s,'lower chord, buckling of the free flange') then c:=36;
      if StrBegin(s,'lower chord, interaction of in-plane flexural buckling and bending') then c:=37;
      if StrBegin(s,'lower chord, global buckling of built-up CS') then c:=38;
      if StrBegin(s,'joint shear') then c:=39;
      if StrBegin(s,'purlinsupport, axial compression and bending - 1') then c:=40;
      if StrBegin(s,'purlinsupport, axial compression and bending - 2') then c:=41;
      if StrBegin(s,'purlinsupport, axial compression and bending - 3') then c:=42;
      if StrBegin(s,'purlinsupport, axial compression and bending - 4') then c:=43;
      if StrBegin(s,'purlinsupport, axial tension and bending - 1') then c:=44;
      if StrBegin(s,'purlinsupport, axial tension and bending - 2') then c:=45;
      if StrBegin(s,'purlinsupport, axial tension and bending - 3') then c:=46;
      if StrBegin(s,'purlinsupport, axial tension and bending - 4') then c:=47;
      if StrBegin(s,'overhang, axial compression and bending - 1') then c:=48;
      if StrBegin(s,'overhang, axial compression and bending - 2') then c:=49;
      if StrBegin(s,'overhang, axial tension and bending - 1') then c:=50;
      if StrBegin(s,'overhang, axial tension and bending - 2') then c:=51;
      if StrBegin(s,'brace end connection, bearing resistance') then c:=52;
      if StrBegin(s,'brace end connection, shear resistance') then c:=53;
      CheckConv:=c;
   end;


begin
   Assign(f,name);
   ReWrite(f);

   //title
   i:=0;
   j:=0;
   repeat
      i:=i+1;
      if name[i]='\' then j:=i;
   until i=Length(name);
   s:=Copy(name,j+1,i-j);
   j1:=Pos('_detail',s);
   s:=Copy(s,1,j1-1);
   WriteLn(f,'*****************************************************************');
   WriteLn(f,' ***************************************************************');
   WriteLn(f);
   WriteLn(f,' DETAILED DOCUMENTATION FOR THE STRUCTURAL DESIGN OF LINDABTRUSS');
   WriteLn(f);
   for i:=1 to ((64-Length(s)) div 2) do
      Write(f,' ');
   WriteLn(f,s);
   WriteLn(f,'                          ',DateTimeToStr(date));
   WriteLn(f);
   WriteLn(f,' ***************************************************************');
   WriteLn(f,'*****************************************************************');
   WriteLn(f);
   WriteLn(f);


   //input data

   WriteLn(f,'**********');
   WriteLn(f);
   WriteLn(f,'INPUT DATA');
   WriteLn(f);
   WriteLn(f,'**********');
   WriteLn(f);

   with mainf do
      begin
         WriteLn(f,'*** Main dimensions');
         WriteLn(f,'  Span:                          ',FloatToStr(EditL.value),' mm');
         if not r_calceave.checked then WriteLn(f,'  Eave height:                   ',FloatToStr(EditVall.value),' mm')
                           else WriteLn(f,'  Eave height:                   ','calculated');
         if not r_calcridge.checked then WriteLn(f,'  Ridge height:                  ',FloatToStr(EditTarej.value),' mm')
                                    else WriteLn(f,'  Ridge height:                  ','calculated');
         if not r_calcSlope.checked
            then begin
                    if BtnSlopePerc.Down then WriteLn(f,'  Slope:                         ',FloatToStr(Round(TanD(EditSlopeDeg.value)*10000)/100),' %');
                    if BtnSlopeDeg.Down then WriteLn(f,'  Slope:                         ',FloatToStr(EditSlopeDeg.value),' deg');
                 end
            else WriteLn(f,'  Slope:                         ','calculated');
         if c_no_eaves.checked then WriteLn(f,'  Eave type:                     ','no eaves');
         if c_stag_bar.checked then WriteLn(f,'  Eave type and size:            ','eave console:  ',FloatToStr(EditEave.value),' mm');
         if c_overhang.checked then WriteLn(f,'  Eave type and size:            ','overhang:  ',FloatToStr(EditEave.value),' mm');
         if c_gorgo_jobb.Checked then WriteLn(f,'  Supports:                      ','hinge on the left, roller on the right');
         if c_gorgo_bal.Checked then WriteLn(f,'  Supports:                       ','roller on the left, hinge on the right');
         if c_nincs_gorgo.Checked then WriteLn(f,'  Supports:                     ','hinges at both supports');
         WriteLn(f,'  Purlin distribution:           ',FloatToStr(EditDist1.value),' mm + ',FloatToStr(EditN.value),' ? ',FloatToStr(EditDist.value),' mm + ',FloatToStr(EditDist2.value),' mm');
         if c_Horizontal.Checked then WriteLn(f,'  Purlin positions are measured: ','horizontally');
         if c_on_roof.Checked then WriteLn(f,'  Purlin positions are measured: ','on the roof');
         WriteLn(f);

         WriteLn(f,'*** Chords, brace members, bolts');
         WriteLn(f,'  Basic upper chord profile:                           ',c_def_fov.text);
         WriteLn(f,'  Upper chord bolt distances (from outside to inside): ',FloatToStr(SpinTopB1.value),' mm + ',FloatToStr(SpinTopRows.value-1),' ? ',FloatToStr(EditTopDist.value),' mm + ',FloatToStr(SpinTopB2.value),' mm');
         if c_upper_purlins.checked then
            begin
               WriteLn(f,'  Lateral supports along upper chord:                  ','at purlins');
               if chkBak.checked then WriteLn(f,'  Purlin support size:                                 ',FloatToStr(EditBak.value),' mm');
            end;
         if c_upper_cont.checked then
            begin
               WriteLn(f,'  Lateral supports along upper chord:                  ','continuous');
               if chkBak.checked then WriteLn(f,'  Load distance at upper chord:                        ',FloatToStr(EditDistUpper.value),' mm');
            end;

         WriteLn(f,'  Basic lower chord profile:                           ',c_def_aov.text);
         WriteLn(f,'  Lower chord bolt distances (from outside to inside): ',FloatToStr(SpinBottomB1.value),' mm + ',FloatToStr(SpinBottomRows.value-1),' ? ',FloatToStr(EditBottomDist.value),' mm + ',FloatToStr(SpinBottomB2.value),' mm');
         if c_lower_cont.checked then WriteLn(f,'  Lateral supports along lower chord:                  ','continuous');
         if c_lower_each.checked then WriteLn(f,'  Lateral supports along lower chord:                  ','each node');
         if c_lower_none.checked then WriteLn(f,'  Lateral supports along lower chord:                  ','none');
         if c_lower_custom.checked then WriteLn(f,'  Lateral supports along lower chord:                ','user defined');
         WriteLn(f,'  Load distance at lower chord:                        ',FloatToStr(EditDistLower.value),' mm');

         WriteLn(f,'  Basic brace member profile:                          ',c_def_racs.text);
         if chkMer.checked then Write(f,'  Columns:                                             perpendicular to upper chord, ');
         if chkFugg.checked then Write(f,'  Columns:                                             vertical, ');
         if StrToInt(EditEvery.text)=1 then Writeln(f,'at every purlin');
         if StrToInt(EditEvery.text)=2 then Writeln(f,'at every 2nd purlin');
         if StrToInt(EditEvery.text)=3 then Writeln(f,'at every3rd purlin');
         if StrToInt(EditEvery.text)=4 then Writeln(f,'at every 4th purlin');
         Write(f,'  Column shift:                                        ');
         if StrToInt(EditSkip.text)=0 then WriteLn(f,'no shift');
         if StrToInt(EditSkip.text)=1 then WriteLn(f,'shift by ',EditSkip.text,' purlin');
         if StrToInt(EditSkip.text)>1 then WriteLn(f,'shift by ',EditSkip.text,' purlins');

         WriteLn(f,'  Bolt diamater and grade:                             ',c_bolts_dia.text,', ',c_bolts_grade.text);

         WriteLn(f);
         WriteLn(f,'*** Load cases');
         for ilc:=1 to nlc do
            begin
               WriteLn(f,'  Load case: ',geom.aLoad[ilc-1].LC,' in load group: ',c_groupmod.items[ c_GroupMod.items.indexof(  geom.LCGroups[ilc-1] )] );

               for irow:=0 to geom.aLoad[ilc-1].lines.count-1 do
                  begin
                     s:=geom.aLoad[ilc-1].lines[irow];
                     StringIn3(s,v);
                     case StrToInt(v[1]) of
                        1: WriteLn(f,'    type: 1, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m');
                        2: WriteLn(f,'    type: 2, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m');
                        3: WriteLn(f,'    type: 3, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m');
                        4: WriteLn(f,'    type: 4, ','acting at: ',v[2]:20,',  p1: ',StrToFloat(v[5]):5:2,' kN/m',', angle: ',StrToFloat(v[9]):5:2,' deg');
                        8: WriteLn(f,'    type: 8, ','acting at: ',v[2]:20,',  F: ',StrToFloat(v[3]):5:1,' kN/m',', angle: ',StrToFloat(v[9]):5:2,' deg',', pos: ',v[10]);
                     end; //case
                  end;
            end;

         WriteLn(f);
         WriteLn(f,'*** Load combinations');
         WriteLn(f,'  Combinations are in accordance with: ',c_env.text);
         WriteLn(f,'  The applied partial and combination factors as follows:');
         for j:=1 to groups.ColCount do
            begin
               Write(f,'  ',STrunc(Groups.cells[0,j],30):30,STrunc(Groups.cells[1,j],15):15);
               if groups.Cells[1,j]=get('Permanent') then
                  begin
                     Write(f,'  ',' gam_G,inf=',Groups.cells[2,j]);
                     Write(f,'  ',' zeta =',Groups.cells[3,j]);
                     Write(f,'  ',' gam_G,sup=',Groups.cells[4,j]);
                     Write(f,'  ',' gam_G,A=',Groups.cells[5,j]);
                //     for i:=2 to 5 do
                //        Write(f,'  ',Groups.cells[i,j]);
                     WriteLn(f);
                  end;
               if groups.Cells[1,j]=get('Variable') then
                  begin
                     Write(f,'  ',' gam_Q    =',Groups.cells[7,j]);
                     Write(f,'  ',' psi_0=',Groups.cells[8,j]);
                     Write(f,'  ',' psi_1    =',Groups.cells[9,j]);
                     Write(f,'  ',' psi_2  =',Groups.cells[10,j]);
                //     for i:=7 to 10 do
                //        Write(f,'  ',Groups.cells[i,j]);
                     WriteLn(f);
                  end;
               if groups.Cells[1,j]=get('Accidental') then
                  begin
                   //  WriteLn(f,'  ',Groups.cells[6,j]);
                     WriteLn(f,'  ',' gam_A=    ',Groups.cells[6,j]);
                  end;
            end;
         WriteLn(f);


         WriteLn(f);
         WriteLn(f,'*** Checkpoints');
         //WriteLn(f,'  Checkpoint density used in calculation: ',AnsiReplaceStr( TMenuitem( Average1 ).Caption, '&',''));
         Write(f,'  Checkpoint density used in calculation: ');
         if StrBegin(BtnCheckPoints.Caption,Get('Checkpoints: Rare')) then WriteLn(f,'rare');
         if BtnCheckPoints.Caption=Get('Checkpoints: Average') then WriteLn(f,'average');
         if BtnCheckPoints.Caption=Get('Checkpoints: Dense') then WriteLn(f,'dense');

      end; //with

   WriteLn(f);
   WriteLn(f);
   WriteLn(f);


   //geometry of static model

   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'NODAL COORDINATES, SUPPORTS');
   WriteLn(f);
   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,' node # ','  x [mm] ','   y [mm] ','  Rx [N/mm] ','  Ry [N/mm] ');
   WriteLn(f);

   for inode:=1 to nnode do
      begin
         Write(f,inode:6,mainf.geom.aStatCsp[inode-1].x:10:2,mainf.geom.aStatCsp[inode-1].y:10:2);
         if mainf.geom.aStatCsp[inode-1].CsuklosTam then WriteLn(f,FixSpring:12:0,FixSpring:12:0);
         if mainf.geom.aStatCsp[inode-1].GorgosTam then WriteLn(f,'           0',FixSpring:12:0);
         if (not mainf.geom.aStatCsp[inode-1].CsuklosTam) and (not mainf.geom.aStatCsp[inode-1].GorgosTam) then WriteLn(f,'          --','          --');
      end;
   WriteLn(f);
   WriteLn(f);


   WriteLn(f,'*********************************');
   WriteLn(f);
   WriteLn(f,'BEAM ELEMENTS OF THE STATIC MODEL');
   WriteLn(f);
   WriteLn(f,'*********************************');
   WriteLn(f);
   WriteLn(f,'  elem ',' start ','   end ','  structural function           ','   profile       ',' single/double ','   area ','  inertia ');
   WriteLn(f,'    no ','  node ','  node ','                                ','                 ','               ','  [cm2] ','    [cm4] ');
   WriteLn(f);

   for ib:=1 to nbeam do
      begin
         Write(f,ib:6);
         inode:=0;

         //start/end node
         Write(f,beamtop[ib,1]:7,beamtop[ib,2]:7);

         //structural function
         j1:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
         case mainf.geom.ARud[j1-1].tipus of
            rAlsoOv:      Write(f,'   lower chord                   ');
            rFelsoOv:     Write(f,'   upper chord                   ');
            rRacs:        Write(f,'   brace member                  ');
            rKulsoOszlop: Write(f,'   brace member (external column)');
            rBelsoOszlop: Write(f,'   brace member (internal column)');
            rBak:         Write(f,'   purlin support                ');
         end;


         //checking if the element is a structural element or not
         nonstruct:=False;
         if mainf.geom.aStatVaz[ib-1].Stagbar then nonstruct:=True;
         if mainf.geom.aStatVaz[ib-1].PurlinSupport and mainf.chkFugg.Checked then nonstruct:=True;
         if mainf.geom.aStatVaz[ib-1].PurlinSupport and mainf.chkMer.Checked then
            begin
               braceconnect:=False;
               for jb:=1 to nbeam do
                  if jb<>ib then
                     begin
                        if (beamtop[jb,1]=beamtop[ib,1]) or (beamtop[jb,1]=beamtop[ib,2])
                           or (beamtop[jb,2]=beamtop[ib,1]) or (beamtop[jb,2]=beamtop[ib,2]) then
                              begin
                                 j1:=mainf.geom.aStatvaz[jb-1].RudIdx+1;
                                 if mainf.geom.ARud[j1-1].tipus in [rBelsoOszlop] then
                                    braceconnect:=True;
                              end;
                     end;
               if not braceconnect then nonstruct:=True;
            end;

         //profile name
         if nonstruct
            then
               begin
                  Write(f,'   ','not defined');
               end
            else
               begin
                  i:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                  Write(f,'   ',mainf.geom.aRud[i-1].profile.profile);
               end;

         //simple/double
         if nonstruct
            then
               begin
                  Write(f,'    ','--');
               end
            else
               begin
                  i:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                  if (mainf.geom.ARud[i-1].tipus in [rRacs,rBelsoOszlop]) and (mainf.geom.ARud[i-1].dupla)
                     then Write(f,'      ','double')
                     else Write(f,'      ','single');
               end;

         //A, I
         if nonstruct
            then Write(f,'         ',RigidA/100:10:0,RigidI/10000:10:0)
            else Write(f,'     ',Av[ib]/100:10:2,Iv[ib]/10000:10:2);

         WriteLn(f);
      end;


   WriteLn(f);
   WriteLn(f);


   //load combinations

   WriteLn(f,'************************');
   WriteLn(f);
   WriteLn(f,'ACTIVE LOAD COMBINATIONS');
   WriteLn(f);
   WriteLn(f,'************************');
   WriteLn(f);

   WriteLn(f);
   WriteLn(f,'*** ULS combinations');
   WriteLn(f);
   //WriteLn(f,'name   ','      permanent load',' + variable loads ');

   for ilc:=1 to mainf.comb.rowcount-1 do
      if (mainf.comb.cells[mainf.comb.colcount-1,ilc]='On') and
         (mainf.comb.cells[mainf.comb.colcount-2,ilc]='ULS') then
         begin
            Write(f,'  ',STrunc(mainf.comb.cells[0,ilc],10):10);
            for i:=1 to mainf.comb.colcount-3 do
               begin
                  if i=1 then Write(f,'   ',mainf.comb.cells[i,ilc])
                         else begin
                                 if (i mod 2)=1 then Write(f,' + ',mainf.comb.cells[i,ilc])
                                                else Write(f,'?','"',STrunc(mainf.comb.cells[i,ilc],10):10,'"');
                              end;
               end;
            WriteLn(f);
         end;

   WriteLn(f);
   WriteLn(f,'*** SLS combinations');
   WriteLn(f);
   //WriteLn(f,'name   ','      permanent load',' + variable loads ');

   for ilc:=1 to mainf.comb.rowcount-1 do
      if (mainf.comb.cells[mainf.comb.colcount-1,ilc]='On') and
         (mainf.comb.cells[mainf.comb.colcount-2,ilc]='SLS') then
         begin
            Write(f,'  ',STrunc(mainf.comb.cells[0,ilc],10):10);
            for i:=1 to mainf.comb.colcount-3 do
               begin
                  if i=1 then Write(f,'   ',mainf.comb.cells[i,ilc])
                         else begin
                                 if (i mod 2)=1 then Write(f,' + ',mainf.comb.cells[i,ilc])
                                                else Write(f,'?','"',STrunc(mainf.comb.cells[i,ilc],10):10,'"');
                              end;
               end;
            WriteLn(f);
         end;

   WriteLn(f);
   WriteLn(f);


   //internal forces, displacements

   WriteLn(f,'**********************************************');
   WriteLn(f);
   WriteLn(f,'INTERNAL FORCES, DISPLACEMENTS FROM LOAD CASES');
   WriteLn(f);
   WriteLn(f,'**********************************************');
   WriteLn(f);
   WriteLn(f,' elem #','   lc # ','  pos [%] ','   N [kN] ','   V [kN] ','  M [kNm] ','  dx [mm] ','  dy [mm] ');
   WriteLn(f);

   for ib:=1 to nbeam do
      begin
         for ilc:=1 to nlc do
            begin
               for icp:=1 to npoints[ib] do
                  begin
                     dx:=defx[ib,ilc,icp];
                     dy:=defy[ib,ilc,icp];
                     LtoG(dx,dy,alfv[ib]);
                     WriteLn(f,ib:7,ilc:7,(100*xrel[ib,icp]):10:1,Nel[ib,ilc,icp]/1000:10:3,Tel[ib,ilc,icp]/1000:10:3,Mel[ib,ilc,icp]/1000000:10:3,dx:10:3,dy:10:3);
                  end; //checkpoint
            end; //beam
         WriteLn(f);
      end;  //load case

   WriteLn(f);
   WriteLn(f);


   //verification results

   for ib:=1 to nbeam do
      begin

         WriteLn(f,'***************************************');
         WriteLn(f);
         WriteLn(f,'CHECK OF ELEMENT#',IntToStr(ib),' OF THE STATIC MODEL');
         WriteLn(f);
         WriteLn(f,'***************************************');
         WriteLn(f);

         umax:=0;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               for icp:=activecp[ib,1] to activecp[ib,2] do
                  for i:=1 to 20 do
                     umax:=MaxReal(umax,ell[ib,ilc,icp,i].util);

         if umax<0.00001
            then
               begin
                  WriteLn(f,'  This element is a non-strucutral element and not checked.');
                  WriteLn(f);
                  WriteLn(f);
               end
            else
               begin
                  j1:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                  WriteLn(f);
                  Write(f,'  Structural function: ');
                  case mainf.geom.ARud[j1-1].tipus of
                     rAlsoOv:      WriteLn(f,'lower chord');
                     rFelsoOv:     WriteLn(f,'upper chord');
                     rRacs:        WriteLn(f,'brace member');
                     rKulsoOszlop: WriteLn(f,'brace member (external column)');
                     rBelsoOszlop: WriteLn(f,'brace member (internal column)');
                     rBak:         WriteLn(f,'purlin support');
                  end;   ;
                  WriteLn(f);
                  yet:=False;
                  for j:=1 to 51 do
                     begin
                        yet2:=False;
                        for ilc:=KombKezd to KombVeg do
                           if aLCOrComb[ilc].uls then
                              for i:=1 to 20 do
                                 begin
                                    if CheckConv(ell[ib,ilc,activecp[ib,1],i].tonkr) = j then
                                       begin
                                          if (not yet) then
                                             begin
                                                WriteLn(f);
                                                WriteLn(f,'*** Checks along the member length');
                                                WriteLn(f);
                                                Write(f,'  Check points (in % of element length):                            ');
                                                for icp:=activecp[ib,1] to activecp[ib,2] do
                                                   Write(f,xrel[ib,icp]*100:6:1);
                                                WriteLn(f);
                                                yet:=True;
                                             end;
                                          if (not yet2) then
                                             begin
                                                WriteLn(f);
                                                WriteLn(f,'  Check type: ',ell[ib,ilc,activecp[ib,1],i].tonkr);
                                                yet2:=True;
                                             end;

                                          Write(f,'    Utilizations in checkpoints (in %), in combination "',STrunc(aLCorComb[ilc].name,10):10,'":');
                                          for icp:=activecp[ib,1] to activecp[ib,2] do
                                             Write(f,ell[ib,ilc,icp,i].util*100:6:1);
                                          WriteLn(f);
                                       end;
                              end; //load combinations
                     end; //check type
                  yet:=False;   
                  for j:=52 to 53 do
                     begin
                        yet2:=False;
                        for ilc:=KombKezd to KombVeg do
                           if aLCOrComb[ilc].uls then
                              for i:=1 to 20 do
                                 begin
                                    if CheckConv(ell[ib,ilc,1,i].tonkr) = j then
                                       begin
                                          if (not yet) then
                                             begin
                                                WriteLn(f);
                                                WriteLn(f,'*** Connection checks');
                                                WriteLn(f);
                                                Write(f,'  Check points (in % of element length):                            ');
                                                WriteLn(f,xrel[ib,1]*100:6:1,xrel[ib,npoints[ib]]*100:6:1);
                                                yet:=True;
                                             end;
                                          if (not yet2) then
                                             begin
                                                WriteLn(f);
                                                WriteLn(f,'  Check type: ',ell[ib,ilc,1,i].tonkr);
                                                WriteLn(f);
                                                yet2:=True;
                                             end;

                                          Write(f,'    Utilizations in checkpoints (in %), in combination "',STrunc(aLCorComb[ilc].name,10):10,'":');
                                          WriteLn(f,ell[ib,ilc,1,i].util*100:6:1,ell[ib,ilc,npoints[ib],i].util*100:6:1);
                                       end;
                                 end; //load combinations
                     end; //check type

                  WriteLn(f);
                  WriteLn(f);

               end; //if structural member

      end;




   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'MAXIMAL UTILIZATIONS IN ULS');
   WriteLn(f);
   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'elem #','   utilization','     combination','   position','    governing check');
   WriteLn(f);

   ULSmax:=-1;
   for ib:=1 to nbeam do
      begin
         umax:=ell[ib,1,1,1].util;
         smax:=ell[ib,1,1,1].tonkr;
         cmax:=KombKezd;
         xmax:=xrel[ib,1];
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  for icp:=1 to npoints[ib] do
                     begin
                        for i:=1 to length(ell[ib,ilc,icp]) do
                        if ell[ib,ilc,icp,i].util>umax then
                           begin
                              umax:=ell[ib,ilc,icp,i].util;
                              smax:=ell[ib,ilc,icp,i].tonkr;
                              cmax:=ilc;
                              xmax:=xrel[ib,icp];
                           end;

                     end; //checkpoint
               end; //load combination
               if umax>ULSmax then ULSmax:=umax;

               if Abs(umax)<0.00001
                  then WriteLn(f,ib:5,'         --    ','    --             ','    --   ','   not checked')
                  else WriteLn(f,ib:5,' ',umax*100:10:2,' %      ',STrunc(aLCorComb[cmax].name,10):10,' ',xmax*100:10:2,' %    ',smax)

         end;  //beam

   WriteLn(f);
   WriteLn(f,'Maximal ULS utilization:  ',ULSmax*100:10:2,' %');

   WriteLn(f);
   WriteLn(f);

   WriteLn(f,'***************************');
   WriteLn(f);
   WriteLn(f,'MAXIMAL UTILIZATIONS IN SLS');
   WriteLn(f);
   WriteLn(f,'***************************');
   WriteLn(f);

   WriteLn(f,'Displacement limit:  L/',mainf.DeflLimit.text);
   WriteLn(f);

   SLSMax:=0;
   for ilc:=KombKezd to KombVeg do
      if aLCOrComb[ilc].sls then
         begin
            if aLCorComb[ilc].MaxSLSChkRatio>SLSMax then SLSMax:=aLCorComb[ilc].MaxSLSChkRatio;
            WriteLn(f,STrunc(aLCorComb[ilc].name,10):10,': ',aLCorComb[ilc].MaxSLSChkRatio:10:2,' %');
         end;
   WriteLn(f);
   WriteLn(f,'Maximal SLS utilization:  ',SLSmax:10:2,' %');
   WriteLn(f);
   WriteLn(f);


   //output for extra checks

   WriteLn(f,'******************************');
   WriteLn(f);
   WriteLn(f,'OUTPUT DATA FOR FURTHER CHECKS');
   WriteLn(f);
   WriteLn(f,'******************************');
   WriteLn(f);

   with mainf do
      begin
         WriteLn(f,'*** Max/min reaction forces');
         for ib:=1 to 2 do
            begin
               Rxmax:=0;
               Rxmin:=0;
               Rymax:=0;
               Rymin:=0;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if reac[ilc,ib,1]>Rxmax then Rxmax:=reac[ilc,ib,1];
                        if reac[ilc,ib,1]<Rxmin then Rxmin:=reac[ilc,ib,1];
                        if reac[ilc,ib,2]>Rymax then Rymax:=reac[ilc,ib,2];
                        if reac[ilc,ib,2]<Rymin then Rymin:=reac[ilc,ib,2];
                     end;
               if ib=1 then Write(f,'  Left support,  ')
                       else Write(f,'  Right support, ');
               WriteLn(f,'horizontal reactions:   ','Rxmin=',Rxmin/1000:8:1,' kN    ','Rxmax=',Rxmax/1000:8:1,' kN');
               if ib=1 then Write(f,'  Left support,  ')
                       else Write(f,'  Right support, ');
               WriteLn(f,'vertical reactions:     ','Rymin=',Rymin/1000:8:1,' kN    ','Rymax=',Rymax/1000:8:1,' kN');
            end;
         WriteLn(f);


         WriteLn(f,'*** Internal forces on lower chord connection');

         //finding lower chord chord connection
         icp:=-1;
         inode:=0;
         repeat
            inode:=inode+1;
            ct:=0;
            for ib:=1 to nbeam do
               if (beamtop[ib,1]=inode) or (beamtop[ib,2]=inode) then
                  begin
                     ct:=ct+1;
                     ib2[ct]:=ib;
                  end;
            if (ct=2) and (mainf.geom.aStatVaz[ib2[1]-1].oldal<>mainf.geom.aStatVaz[ib2[2]-1].oldal) then
               begin
                  ib:=ib2[1];
                  j1:=mainf.geom.aStatvaz[ib2[1]-1].RudIdx+1;
                  j2:=mainf.geom.aStatvaz[ib2[2]-1].RudIdx+1;
                  if (mainf.geom.ARud[j1-1].tipus in [rAlsoOv]) and (mainf.geom.ARud[j2-1].tipus in [rAlsoOv]) then
                     begin
                        if beamtop[ib,1]=inode then icp:=1;
                        if beamtop[ib,2]=inode then icp:=npoints[ib];
                      end;
               end;
         until (icp>0) or (i=nnode);

         //internal forces: node=inode, beam=ib, checkpoint=icp

         //max N
         Nmax:=-1000000000;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Nel[ib,ilc,icp]>Nmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Max N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
         //min N
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Nel[ib,ilc,icp]<Nmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Min N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

         //max T
         Tmax:=-1000000000;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Tel[ib,ilc,icp]>Tmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Max V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
         //min T
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Tel[ib,ilc,icp]<Tmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Min V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

         //max M
         Mmax:=-1000000000;
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Mel[ib,ilc,icp]>Mmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Max M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
         //min M
         for ilc:=KombKezd to KombVeg do
            if aLCOrComb[ilc].uls then
               begin
                  if Mel[ib,ilc,icp]<Mmax then
                      begin
                         Nmax:=Nel[ib,ilc,icp];
                         Tmax:=Tel[ib,ilc,icp];
                         Mmax:=Mel[ib,ilc,icp];
                      end;
               end; //load comb
         WriteLn(f,'  Min M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

         WriteLn(f);


         WriteLn(f,'*** Internal forces on ridge connection');

         //finding ridge
         icp:=-1;
         inode:=0;
         repeat
            inode:=inode+1;
            ct:=0;
            for ib:=1 to nbeam do
               if (beamtop[ib,1]=inode) or (beamtop[ib,2]=inode) then
                  begin
                     ct:=ct+1;
                     ib2[ct]:=ib;
                  end;
            if (ct=2) and (mainf.geom.aStatVaz[ib2[1]-1].oldal<>mainf.geom.aStatVaz[ib2[2]-1].oldal) then
               begin
                  j1:=mainf.geom.aStatvaz[ib2[1]-1].RudIdx+1;
                  j2:=mainf.geom.aStatvaz[ib2[2]-1].RudIdx+1;
                  if (mainf.geom.ARud[j1-1].tipus in [rFelsoOv]) and (mainf.geom.ARud[j2-1].tipus in [rFelsoOv]) then icp:=1;
               end;
         until (icp>0) or (i=nnode);


         for i:=1 to 2 do
            begin
               ib:=ib2[i];
               if beamtop[ib,1]=inode then icp:=1;
               if beamtop[ib,2]=inode then icp:=npoints[ib];

               //internal forces: node=inode, beam=ib, checkpoint=icp

               if mainf.geom.aStatVaz[ib-1].oldal=oBal then WriteLn(f,'  Internal forces from the left');
               if mainf.geom.aStatVaz[ib-1].oldal=oJobb then WriteLn(f,'  Internal forces from the right');

               //max N
               Nmax:=-1000000000;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Nel[ib,ilc,icp]>Nmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Max N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
               //min N
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Nel[ib,ilc,icp]<Nmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Min N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

               //max T
               Tmax:=-1000000000;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Tel[ib,ilc,icp]>Tmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Max V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
               //min T
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Tel[ib,ilc,icp]<Tmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Min V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

               //max M
               Mmax:=-1000000000;
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Mel[ib,ilc,icp]>Mmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Max M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
               //min M
               for ilc:=KombKezd to KombVeg do
                  if aLCOrComb[ilc].uls then
                     begin
                        if Mel[ib,ilc,icp]<Mmax then
                            begin
                               Nmax:=Nel[ib,ilc,icp];
                               Tmax:=Tel[ib,ilc,icp];
                               Mmax:=Mel[ib,ilc,icp];
                            end;
                     end; //load comb
               WriteLn(f,'    Min M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

               WriteLn(f);
            end; //for i:=1 to 2

            
         WriteLn(f,'*** Internal forces on eaves');

         if c_no_eaves.checked
            then
               begin
                  WriteLn(f,'  There are no eaves.');
               end

         else
            begin
               for i:=1 to 2 do
                  begin
                     case i of
                        1: side:=oBal;
                        2: side:=oJobb;
                     end;

                     //finding left/right eave
                     if c_overhang.checked then     //overhang
                        begin
                           ib:=0;
                           repeat
                              ib:=ib+1;
                              j1:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                           until (geom.ARud[j1-1].tipus=rFelsoov) and (geom.aStatVaz[ib-1].oldal=side)
                              and (geom.aStatVaz[ib-1].Overhang) and (geom.aStatVaz[ib-1].tipus=rFelsoOv)
                              and (CheckParallel(geom.ARud[j1-1].l,geom.aStatVaz[ib-1].v)) or (ib=nbeam);
                        end; //overhang
                     if c_stag_bar.checked then     //stag bar
                        begin
                           ib:=0;
                           repeat
                              ib:=ib+1;
                              j1:=mainf.geom.aStatvaz[ib-1].RudIdx+1;
                           until (geom.ARud[j1-1].tipus=rFelsoov) and (geom.aStatVaz[ib-1].oldal=side)
                              and (geom.aStatVaz[ib-1].Stagbar) and (geom.aStatVaz[ib-1].tipus=rFelsoOv)
                              and (CheckParallel(geom.ARud[j1-1].l,geom.aStatVaz[ib-1].v)) or (ib=nbeam);
                        end; //overhang)

                     case i of
                        1: j2:=1;
                        2: j2:=-1;
                     end;
                     if (j2*mainf.geom.aStatVaz[ib-1].v.x1>j2*mainf.geom.aStatVaz[ib-1].v.x2)
                        then icp:=1
                        else icp:=npoints[ib];


                     //internal forces: beam=ib, checkpoint=icp

                     case i of
                        1: WriteLn(f,'  Internal forces on the the left eave');
                        2: WriteLn(f,'  Internal forces on the the right eave');
                     end;

                     //max N
                     Nmax:=-1000000000;
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Nel[ib,ilc,icp]>Nmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Max N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
                     //min N
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Nel[ib,ilc,icp]<Nmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Min N, simultaneous V and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

                     //max T
                     Tmax:=-1000000000;
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Tel[ib,ilc,icp]>Tmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Max V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
                     //min T
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Tel[ib,ilc,icp]<Tmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Min V, simultaneous N and M:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

                     //max M
                     Mmax:=-1000000000;
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Mel[ib,ilc,icp]>Mmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Max M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');
                     //min M
                     for ilc:=KombKezd to KombVeg do
                        if aLCOrComb[ilc].uls then
                           begin
                              if Mel[ib,ilc,icp]<Mmax then
                                  begin
                                     Nmax:=Nel[ib,ilc,icp];
                                     Tmax:=Tel[ib,ilc,icp];
                                     Mmax:=Mel[ib,ilc,icp];
                                  end;
                           end; //load comb
                     WriteLn(f,'    Min M, simultaneous N and V:   ',' N=',Nmax/1000:7:2,' kN, ',' V=',Tmax/1000:7:2,' kN, ',' M=',Mmax/1000000:7:2,' kNm');

                     WriteLn(f);
                  end; //for i:=1 to 2

            end; //if there are eaves

      end;   //with







   Close(f);

end; //StatDocFull




procedure CalcStructure(var beamtop:beamtopoltype;var Av,Iv:realbeamtype;
                        var npoints:intbeamtype;var xrel:realbeampointstype;
                        var Nel,Tel,Mel,defx,defy:realbeamlcpointstype;
                        var reac:reallcnodeproptype;var activecp:beamtopoltype);

//this routine calculates the deflections and inernal forces for the whole structure
//  for all the load cases

//output
//
//npoints             vector of integers, its i-th element gives how many checkpoints
//                        are defined within the i-th element
//                        (element numbering: according to statical model)
//xrel                2D matrix
//                       1st dimension: number of beam element
//                       2nd dimension: positions of the checkpoints along the beam element
//                            position is given as relative to beam element length, i.e., position is between 0 and 1
//                            (can be equal to 0 or 1, too)
//Nel,Tel, Mel        3D matrices
//                       1st dimension: number of beam element
//                       2nd dimension: number of load case
//                       3rd dimension: the value of internal force (N, T or M) in the checkpoints
//                          (e.g., Nel[1,2,3] gives the
//                                  - normal force
//                                  - at the 3rd checkpoint
//                                  - of the 1st beam element
//                                  - in the 2nd load case )
//defx,defy           3D matrices
//                       1st dimension: number of beam element
//                       2nd dimension: number of load case
//                       3rd dimension: the value of deflection (in the x or y direction) in the checkpoints
//                          Note: x and y are in local system, x: parallel, y: perpendicular
//                          (e.g., defy[4,5,6] gives the
//                                  - displacement perpendicular to the beam
//                                  - at the 6th checkpoint
//                                  - of the 4th beam element
//                                  - in the 5th load case )
//reac                3D matrix
//                       1st dimension: number of load case
//                       2nd dimension: number of support (node)
//                       3rd dimension: the value of reactions (1:x dir force, 2:y dir force, 3:moment)
//                          Note: x and y are in global system, x: horizontal, y: vertical
//                          (e.g., reac[3,2,1] gives the
//                                  - horizontal reaction force
//                                  - at the 2nd support
//                                  - in the 3rd load case



var
   nlc,nbeam,nnode,nsup,nwidth,i,j,k:  integer;
   E,alft,Lmin:   real;
   Lv,alfv:    realbeamtype;
   Kglob:  stifmattype;
   nodes,supports: nodeproptype;
   nloads:    intlcbeamtype;
   loadtyps:   intlcbeamloadtype;
   loadprops:   multiloadproptype;
   qglob,xglob:   reallcdoftype;
   intforces:     real3type;
   defs:    real2type;
   supploads:reallcnodeproptype;
   f:text;
   offset:realbeam2type;
   tonew,toold: intnodetype;

begin
   ReadInpData(nbeam,nnode,nsup,nodes,beamtop,supports,supploads,E,alft,Av,Iv,nlc,nloads,loadtyps,loadprops,tonew,toold);
   //Ez sz?ks?ges a val?di terhek ?tvitel?hez
   mainf.geom.TransferLoads(nlc,nloads,loadtyps,loadprops,supploads);
   //Ez sz?ks?ges a teherinf? megjelen?t?s?hez
   mainf.geom.TransferLoads(nlc_,nloads_,loadtyps_,loadprops_,supploads_);

   nwidth:=BandCalc(nbeam,beamtop);
   BeamProps(nnode,nbeam,nodes,beamtop,Lv,alfv);
   Lmin:=SetLMin(nbeam,Lv);
   BeamEndOffset(nbeam,beamtop,Lv,alfv,offset);
   Kglobal(nbeam,nnode,nwidth,nsup,beamtop,supports,Av,Iv,Lv,alfv,E,Kglob);
   Calcqglob(nlc,nsup,nloads,loadtyps,loadprops,nbeam,beamtop,Lv,Av,alfv,supports,supploads,E,alft,qglob);
   EquaSolve(nnode*3,nwidth,nlc,Kglob,qglob,xglob);
   CalcReactions(nloads,loadtyps,loadprops,nlc,nsup,supports,supploads,xglob,reac);


   //Assign(f,'igv_bak.txt');
   //ReWrite(f);

   for i:=1 to nbeam do
      begin
         BeamCheckPoints(nloads,loadtyps,loadprops,nlc,i,Lv[i],Lmin,offset[i],npoints[i],xrel[i]);
         CalcIntForces(nloads,loadtyps,loadprops,nlc,i,beamtop,Lv,Av,Iv,alfv,E,alft,xglob,npoints[i],xrel[i],Nel[i],Tel[i],Mel[i]);
         ReduceCheckPoints(nlc,Lv[i],offset[i],npoints[i],xrel[i],Nel[i],Tel[i],Mel[i],activecp[i]);
         CalcDeformations(nloads,loadtyps,loadprops,nlc,i,beamtop,Lv,Av,Iv,alfv,E,alft,xglob,npoints[i],xrel[i],defx[i],defy[i]);

   //      if mainf.geom.aStatVaz[i-1].PurlinSupport then
   //          begin
   //             Writeln(f,'Elem sorsz?ma a stat. v?zban: ',i:5);
   //             k:=npoints[i];
   //             for j:=1 to nlc do
   //                begin
   //                   Writeln(f,'     Tehereset sorsz?ma: ',j:5);
   //                   Writeln(f,'          Kezd?pont ig?nybev?telei (N,T,M): ',Nel[i,j,1]/1000:15:3,Tel[i,j,1]/1000:15:3,Mel[i,j,1]/1000000:15:3);
   //                   Writeln(f,'          V?gpont ig?nybev?telei   (N,T,M): ',Nel[i,j,k]/1000:15:3,Tel[i,j,k]/1000:15:3,Mel[i,j,k]/1000000:15:3);
   //                end;
   //          end;

      end;

   //reset node numbering
   for i:=1 to nbeam do
      begin
         beamtop[i,1]:=tonew[beamtop[i,1]];
         beamtop[i,2]:=tonew[beamtop[i,2]];
      end;
   //Close(f);

   //CalcOutShort(nlc,nbeam,npoints,xrel,'out_short.txt',Nel,Tel,Mel,defx,defy);
   //CalcOutFull(nlc,nbeam,npoints,xrel,'out_full.txt',Nel,Tel,Mel,defx,defy);

   FormProgress.visible:=false;
   if mainf.geom.ErrList.Count>0 then
     MessageDlg(get('Warning!')+LF+mainf.geom.ErrList.gettext,mtWarning,[mbok],0);
end; //CalcStructure



end.


