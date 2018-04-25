ptsorig = [[-1,0.25],[-1.1,0.15],[-1.2,0],[-1.1,0],[-1,0],[-0.9,0],[-0.8,0],[-0.7,0],[-0.6,0],[-0.5,0.008],[-0.4,0.016],[-0.3,0.024],[-0.2,0.032],[-0.1,0.041],[0,0.05],[0.1,0.075],[0.2,0.1],[0.3,0.125],[0.4,0.14],[0.5,0.16],[0.6,0.2],[0.7,0.25],[0.8,0.3],[0.9,0.35],[1,0.4],[1.1,0.475],[1.2,0.6],[1.225,0.7],[1.2,0.8],[1.1,0.83],[1,0.84],[0.9,0.85],[1,0.95],[1.1,1.05],[1.15,1.1],[1.1,1.1],[1,1.1],[0.9,1.1],[0.8,1.1],[0.7,1.1],[0.6,1.1],[0.5,1.09],[0.4,1.08],[0.3,1.07],[0.2,1.06],[0.1,1.05],[0,1.04],[-0.1,1.03],[-0.2,1.02],[-0.3,1.01],[-0.4,1],[-0.5,0.975],[-0.6,0.95],[-0.7,0.925],[-0.8,0.9],[-0.9,0.865],[-1,0.835],[-1.1,0.8],[-1.2,0.7],[-1.3,0.6],[-1.37,0.5],[-1.37,0.4],[-1.3,0.3],[-1.2,0.28],[-1.1,0.25],[-0.95,0.3],[-1,0.3],[-1.1,0.3],[-1.15,0.3],[-1.1,0.4],[-1,0.4],[-0.9,0.4],[-0.8,0.4],[-0.7,0.4],[-0.6,0.4],[-0.53,0.45],[-0.6,0.5],[-0.7,0.5],[-0.8,0.5],[-0.9,0.53],[-0.95,0.6],[-0.93,0.7],[-0.82,0.8],[-0.7,0.8],[-0.6,0.8],[-0.5,0.8],[-0.4,0.8],[-0.3,0.8],[-0.2,0.8],[-0.25,0.7],[-0.15,0.7],[-0.1,0.75],[0,0.8],[0.1,0.8],[0.2,0.8],[0.3,0.8],[0.4,0.8],[0.5,0.8],[0.55,0.8],[0.5,0.7],[0.4,0.7],[0.3,0.7],[0.2,0.7],[0.1,0.7],[-0.05,0.6],[0.1,0.6],[0.2,0.6],[0.3,0.6],[0.4,0.6],[0.45,0.6],[0.4,0.5],[0.3,0.5],[0.2,0.5],[0.1,0.5],[0,0.5],[-0.1,0.5],[-0.12,0.47],[-0.12,0.43],[-0.07,0.4],[0,0.4],[0.1,0.4],[0.2,0.4],[0.3,0.4],[0.35,0.4],[0.3,0.3],[0.4,0.3],[0.45,0.4],[0.5,0.5],[0.55,0.6],[0.6,0.7],[0.65,0.8],[0.7,0.8],[0.8,0.8],[0.9,0.8],[1,0.8],[1.05,0.8],[1.1,0.77],[1.14,0.7],[1.13,0.6],[1.1,0.55],[1,0.5],[0.9,0.5],[0.8,0.5],[0.7,0.5],[0.75,0.6],[0.8,0.6],[0.9,0.6],[1,0.6],[1.07,0.63],[1.07,0.67],[1,0.7],[0.9,0.7],[0.8,0.7],[0.75,0.7],[0.7,0.6],[0.65,0.5],[0.6,0.4],[0.55,0.3],[0.4,0.3],[0.3,0.3],[0.2,0.3],[0.1,0.3],[0,0.3],[-0.1,0.3],[-0.2,0.36],[-0.23,0.4],[-0.27,0.5],[-0.23,0.6],[-0.2,0.65],[-0.15,0.7],[-0.3,0.7],[-0.4,0.7],[-0.5,0.7],[-0.6,0.7],[-0.7,0.7],[-0.77,0.67],[-0.77,0.63],[-0.7,0.6],[-0.6,0.6],[-0.5,0.6],[-0.43,0.53],[-0.4,0.45],[-0.45,0.35],[-0.5,0.3],[-0.6,0.3],[-0.7,0.3],[-0.8,0.3],[-0.9,0.3],[-0.95,0.3]];
SCALE_FACTOR = 0.05;
function scalepts(i,pts,factor) = (i == -1 ? [] : concat(scalepts(i-1,pts,factor),[[pts[i][0]*factor,pts[i][1]*factor]]));
pts = scalepts(len(ptsorig)-1,ptsorig,SCALE_FACTOR);
d = sqrt(2);
inc = 360/len(pts);
function calcpt1(pt) = sqrt(pow(pt[0]+1,2)+pow(pt[1]+1,2)) - (sqrt(2)-1);
function calcdeg1(pt) = atan((pt[1]+sqrt(2)/2)/(pt[0]+sqrt(2)/2)) - 45;
function generatepts1(i,pts,deg) = (i == len(pts) ? [] : concat(generatepts1(i+1,pts,deg+inc),[[calcpt1(pts[i])*cos(deg+calcdeg1(pts[i])),calcpt1(pts[i])*sin(deg+calcdeg1(pts[i]))]]));
function calcpt2(pt) = sqrt(pow(-pt[0]+1,2)+pow(pt[1]+1,2)) - (sqrt(2)-1);
function calcdeg2(pt) = atan((pt[1]+sqrt(2)/2)/(-pt[0]+sqrt(2)/2)) - 45;
function generatepts2(i,pts,deg) = (i == len(pts) ? [] : concat(generatepts2(i+1,pts,deg-inc),[[calcpt2(pts[i])*cos(deg-calcdeg2(pts[i])),calcpt2(pts[i])*sin(deg-calcdeg2(pts[i]))]]));
c1pts = generatepts1(0,pts,-$t*inc*len(pts));
c2pts = generatepts2(0,pts,$t*inc*len(pts) + 180);
//translate([0,sqrt(2)/2,0]) color("Green") linear_extrude(height=0.001)polygon(points=pts);
module Gear(radius,teeth,teeth_height) {
    for(n=[0:teeth-1]) {
        rotate(a=[0,0,n*360/teeth]){
            linear_extrude(height=0.2) {
                polygon(points=[[0,0],[radius*sin(180/teeth),radius*cos(180/teeth)],[(radius+teeth_height)*sin(135/teeth),(radius+teeth_height)*cos(135/teeth)],[(radius+teeth_height)*sin(45/teeth),(radius+teeth_height)*cos(45/teeth)],[0,radius],[-radius*sin(180/teeth),radius*cos(180/teeth)]]);
            };
        };
    };
}
s=10;
scale(v=[s,s,s]) {
    
    //Disc 1
    difference(){
        translate([-d/2,0,0]) color("Red") linear_extrude(height=0.2) polygon(points=c1pts,convexity=100);
        translate([-d/2,0,-1]) linear_extrude(height=2) rotate([0,0,-$t*360]) polygon(points=[[-0.2,0.1],[0.1,0.1],[0.1,-0.1],[0,-0.1],[0,0],[-0.2,0]]);
    }

    //Disc 2
    difference(){
        translate([d/2,0,0.20]) color("Blue") linear_extrude(height=0.2)polygon(points=c2pts,convexity=100);
        translate([d/2,0,-1]) linear_extrude(height=2) rotate([0,0,$t*360]) polygon(points=[[0.2,0.1],[-0.1,0.1],[-0.1,-0.1],[0,-0.1],[0,0],[0.2,0]]);
    }

    //Gear/axle 1
    union(){
        color("Purple") translate([-d/2,0,0.4]) rotate([0,0,-$t*360]) Gear(0.65,15,0.1);
        translate([-d/2,0,-5]) linear_extrude(height=10) rotate([0,0,-$t*360]) polygon(points=[[-0.19,0.09],[0.09,0.09],[0.09,-0.09],[0.01,-0.09],[0.01,0.01],[-0.19,0.01]]);
        translate([-d/2,0,-5]) linear_extrude(height=1) circle(0.25,$fn=100);
        translate([-d/2,0,4]) linear_extrude(height=1) circle(0.25,$fn=100);
    }

    //Gear/axle 2
    union() {
        color("Yellow") translate([d/2,0,0.4]) rotate([0,0,$t*360]) Gear(0.65,15,0.1);
        translate([d/2,0,-5]) linear_extrude(height=10) rotate([0,0,$t*360]) polygon(points=[[0.19,0.09],[-0.09,0.09],[-0.09,-0.09],[-0.01,-0.09],[-0.01,0.01],[0.19,0.01]]);
        translate([d/2,0,-5]) linear_extrude(height=1) circle(0.25,$fn=100);
        translate([d/2,0,4]) linear_extrude(height=1) circle(0.25,$fn=100);
    }

    //Main scaffolding
    union() {
        difference() {
            translate([-d/2,0,-5.5]) linear_extrude(height=1.5) circle(0.4, $fn=100);
            translate([-d/2,0,-5.05]) linear_extrude(height=1.05) circle(0.26, $fn=100);
        };
        difference() {
            translate([d/2,0,-5.5]) linear_extrude(height=1.5) circle(0.4, $fn=100);
            translate([d/2,0,-5.05]) linear_extrude(height=1.05) circle(0.26, $fn=100);
        };
        difference() {
            translate([-d/2,0,4]) linear_extrude(height=1.5) circle(0.4, $fn=100);
            translate([-d/2,0,4]) linear_extrude(height=1.05) circle(0.26, $fn=100);
            translate([-d/2-0.5,0,3.5]) linear_extrude(height=2) square(1,1);
        };
        difference() {
            translate([d/2,0,4]) linear_extrude(height=1.5) circle(0.4, $fn=100);
            translate([d/2,0,4]) linear_extrude(height=1.05) circle(0.26, $fn=100);
            translate([d/2-0.5,0,3.5]) linear_extrude(height=2) square(1,1);
        };
        translate([0,0,-6.5]) linear_extrude(height=1) difference() { 
            polygon(points=[[-3,0.5],[-3,-3.5],[3,-3.5],[3,0.5]]);
            polygon(points=[[-2,-0.5],[-2,-2.5],[2,-2.5],[2,-0.5]]);
        };
        translate([0,0,5]) linear_extrude(height=1) difference() {
            polygon(points=[[-3,0.5],[-3,-3.5],[3,-3.5],[3,0.5]]);
            polygon(points=[[-2,-0.5],[-2,-2.5],[2,-2.5],[2,-0.5]]);
        };
        
        translate([2,-3.5,-5.5]) linear_extrude(height=10.5) square(1,1);
        translate([-3,-3.5,-5.5]) linear_extrude(height=10.5) square(1,1);
        translate([-0.25,-3.5,-5.5]) linear_extrude(height=10.5) square(0.5,0.5);
        difference(){
            union() {
                translate([-3,-2.5,-0.25]) linear_extrude(height=0.65) square(size=[1,2.5]);
                translate([-2.5,0,-0.25]) linear_extrude(height=0.65) circle(0.5, $fn=100);
            }
            translate([-2.5,0,-0.5]) linear_extrude(height=1) circle(0.2, $fn=100);
        }
        translate([-2.5,0,0.41]) Gear(1.03,23,0.1);
        translate([-2.5,0,-0.5]) linear_extrude(height=1) circle(0.19, $fn=100);
        translate([-2.5,0,-0.3]) linear_extrude(height=0.04) circle(0.3,$fn=100);
        translate([-2.5,0,-0.8]) linear_extrude(height=0.3) circle(0.4,$fn=100);
        difference() {
            union(){
                translate([-2.5,0,-0.8]) linear_extrude(height=0.3) /*rotate([0,0,90]) */polygon(points=[[-0.4,0],[0.4,0],[0.2,-1.4],[-0.2,-1.4]]);
                translate([-2.5,-1.4,-0.8]) linear_extrude(height=0.3) circle(0.2,$fn=100);
            };
            translate([-2.5,-1.375,-1]) linear_extrude(height=1) circle(0.125,$fn=100);
        };
        translate([-2.5,-1.375,-0.49]) linear_extrude(height=0.1) circle(0.2, $fn=100);
        translate([-2.5,-1.375,-0.81]) linear_extrude(height=0.32) circle(0.115, $fn=100);
        translate([-2.5,-1.375,-1.81]) linear_extrude(height=1) circle(0.25, $fn=100);
    }
}


