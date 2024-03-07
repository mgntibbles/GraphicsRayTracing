let sphereCenters = [];
let R = 0;
let sphereAmbients = [];
let sphereDiffuses = [];
let sphereSpeculars = [];
let sphereShiny = [];

let lightCenter = [];
let lightAmb = [];
let pointLight = [];

let numSpheres = 0;

const width = 512;
const height = 512;


function read_file(input){
    let file = input.files[0]
    let reader = new FileReader();
    reader.readAsText(file);

    reader.onload = function() {
        let type = '';
        let contents = this.result.split('\n');
        for (let i =0; i<contents.length; i++){
            let line = contents[i].split(' ');
            type = line[0];
            line = line.filter(function(con){
                return con != '';
            })
            if (type == "light"){
                lightCenter.push(parseFloat(line[1]), parseFloat(line[2]), parseFloat(line[3]));
                lightAmb.push(parseFloat(line[4]), parseFloat(line[5]), parseFloat(line[6]));
                pointLight.push(parseFloat(line[7]), parseFloat(line[8]), parseFloat(line[9]));
            } else if (type == "sphere"){
                sphereCenters.push(parseFloat(line[1]), parseFloat(line[2]), parseFloat(line[3]));
                R = parseFloat(line[4]);
                sphereAmbients.push(parseFloat(line[5]), parseFloat(line[6]), parseFloat(line[7]));
                sphereDiffuses.push(parseFloat(line[8]), parseFloat(line[9]), parseFloat(line[10]));
                sphereSpeculars.push(parseFloat(line[11]), parseFloat(line[12]), parseFloat(line[13]));
                sphereShiny.push(parseFloat(line[14]));
            } else if (type!=''){
                numSpheres = type;
            }
            
        }
        /*
        console.log("light center "+lightCenter);
        console.log("light amb "+lightAmb);
        console.log("point light "+pointLight);
        console.log("sphere centers "+sphereCenters);
        console.log("sphere amb "+sphereAmbients);
        console.log("sphere diffuse "+sphereDiffuses);
        console.log("sphere spec "+sphereSpeculars);
        console.log("sphere shiny "+ sphereShiny);

        console.log(numSpheres);
        */
       calculate_image();
    }
}

function normalize_vector(x, y, z){
    length = Math.sqrt(x**2+y**2+z**2);
    return [x/length, y/length, z/length];
}

function find_discriminant(A, B, C){
    dis = B**2 - 4*A*C;
    return dis;
}

function find_distance(p1, p2){
    distance = Math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2);
    return distance;
}

function dot_product(V1, V2){
    prod = (V1[0]*V2[0]) + (V1[1]*V2[1]) + (V1[2]*V2[2]);
    return prod;
}

function cross_product(V1, V2){
    prod= [(V1[1]*V2[2])-(V1[2]*V2[1]), (V1[2]*V2[0])-(V1[0]*V2[2]), (V1[0]*V2[1])-(V1[1]*V2[0])];
    return prod;
}

function calculate_normal(i, c, Sr, s){
    N = [(i[0] - c[(s*3)])/Sr, (i[1] - c[(s*3)+1])/Sr, (i[2] - c[(s*3)+2])/Sr];
    prd = dot_product(N, [0,0,0]);
    if (prd>=0){
        N = [-N[0], -N[1], -N[2]];
    }
    return normalize_vector(N[0], N[1], N[2]);
}

function calculate_light(loc, i){
    L = [(loc[0]-i[0]), (loc[1]-i[1]), (loc[2]-i[2])];
    return normalize_vector(L[0], L[1], L[2]);
}


function calculate_reflection(N, L){
    con = 2 * (dot_product(N, L));
    temp = [N[0]*con, N[1]*con, N[2]*con];
    Rvec = [temp[0]-L[0], temp[1]-L[1], temp[2]-L[2]];
    return normalize_vector(Rvec[0], Rvec[1], Rvec[2]);
}

function calculate_viewing(i){
    Vv = [0-i[0], 0-i[1], 0-i[2]];
    return normalize_vector(Vv[0], Vv[1], Vv[2]);
}

function illumination_equation(Nv, Lv, Rv, V, sAmb, sDiff, sSpec, lAmb, pLight, s){
    colours = ['r', 'g', 'b'];
    newColour = [];
    let NL = dot_product(Nv, Lv);
    let RV = dot_product(Rv, V);
    if (RV < 0){
        console.log(Nv, Lv);
    }
    for (let i = 0; i < colours.length; i++){
        amb = lAmb[i]*sAmb[i];
        diff = pLight[i]*sDiff[(s*3)+i]*NL;
        spec = pLight[i]*sSpec[(s*3)+i]* (RV)**10.0;
        let col = amb + diff + spec;
        //let col = amb + diff;
        newColour.push(col*255);
        //console.log("amb: "+ amb+" diff: "+diff+" spec: "+spec);
    }
    return newColour;
}

//viewing plane where z = -1
function calculate_image(){
    const canvas = document.querySelector('.myCanvas');
    const ctx = canvas.getContext('2d');

    ctx.fillStyle = 'rgb(0, 0, 0)';
    ctx.fillRect(0, 0, width, height);

    let z = -1
    let vp = [0,0,0];
    let w2 = 2/width;
    let h2 = 2/height;
    let x, y;
    
    for (x=0; x<width-1; x+=1){
        for (y=0; y<height-1;y+=1){
            let ri = []
            //calculate a ray from the viewpoint (0,0,0) for this pixel.
            //Normalize it
            let V = [0-((x*w2)-1), 0-((y*h2)-1),0-z];
            //let V = [((x*w2)-1), ((y*h2)-1), z];
            V = normalize_vector(V[0], V[1], V[2]);
            //calculate if the ray intersects the sphere
            for (s=0; s<numSpheres; s++){
                let A = 1.0
                let B = 2 * ( V[0] * (vp[0] - sphereCenters[(s*3)]) + V[1] * (vp[1] - sphereCenters[(s*3)+1]) + V[2] * (vp[2] - sphereCenters[(s*3)+2]));
                let C  = (vp[0] - sphereCenters[(s*3)])**2 + (vp[1] - sphereCenters[(s*3)+1])**2 + (vp[2] - sphereCenters[(s*3)+2])**2 - R**2;
                dis = find_discriminant(A, B, C);
                
                //console.log("Sphere: "+sphereCenters[(s*3)]+" "+sphereCenters[(s*3)+1]+" "+sphereCenters[(s*3)+2]);
                //console.log("V[0]: "+ V[0]+" vp[0]: "+vp[0]+ " sphereCenters[(s*3)]: "+sphereCenters[(s*3)]+" V[1]: "+ V[1]+ " vp[0]: "+vp[0]+" sphereCenters[(s*3)+1]: "+ sphereCenters[(s*3)+1]+" V[2]: "+V[2]+" vp[2]: "+vp[2]+" sphereCenters[(s*3)+2]: "+sphereCenters[(s*3)+2]+" Sr: "+R);

                //console.log("A: " + A + " B: " +B + " C: "+ C);
                if (dis < 0.01 && dis > -0.01){
                    t0 = (-B - Math.sqrt(dis))/(2*A);
                    ri = [([vp[0]] + V[0]*t0), ([vp[1]] + V[1]*t0), (vp[2] + V[2]*t0)];
                    //console.log("check "+x+" "+y);
                    //console.log(ri.length);
                } else if (dis > 0) {
                    t0 = (-B - Math.sqrt(dis))/(2*A);
                    t1 = (-B + Math.sqrt(dis))/(2*A);
                    let ri0 = [(vp[0] + V[0]*t0), (vp[1] + V[1]*t0), (vp[2] + V[2]*t0)];
                    let ri1 = [(vp[0] + V[0]*t1), (vp[1] + V[1]*t1), (vp[2] + V[2]*t1)];
                    //save closest intersection value
                    distance1 = find_distance(vp, ri0);
                    distance2 = find_distance(vp, ri1);
                    if (distance1>distance2){
                        ri = ri0;
                    } else {
                        ri = ri1;
                    }
                    //ctx.fillStyle = `rgb(100,100,100)`;
                    //ctx.fillRect(x, y, 1, 1);
                    //console.log("check"+x+" "+y);
                    //console.log(ri.length);
                }
                if (ri.length == 3){
                    Vv = calculate_viewing(ri);
                    Nv = calculate_normal(ri, sphereCenters, R, s);
                    Lv = calculate_light(lightCenter, ri);
                    Rv = calculate_reflection(Nv, Lv); //- not needed for this assignment?
                    colour = illumination_equation(Nv, Lv, Rv, Vv, sphereAmbients, sphereDiffuses, sphereSpeculars, lightAmb, pointLight, s);
                    rgb = `rgb(`+colour[0]+`,`+colour[1]+`,`+colour[2]+`)`;
                    //console.log(rgb);
                    ctx.fillStyle =  rgb;
                    ctx.fillRect(x, y, 1, 1);
                } 
                
            }
            //calculate if the ray intersects the object
                //if yes find closest intersection on the object and determine if intersection is closer than the closest sphere intersection
                //        if it was closer than sphere, use .obj
                //calculate the normal if not provided in .obj fil
            //if the ray intersects the sphere or obj
            
                //calculate the reflection vector
                //calculate the light vector
                //calculate the illumination value for the pixel 
                //draw on screen
        }
    }
}

function fill(){
    const canvas = document.querySelector('.myCanvas');

    // set size of 2D image


    //const width = 1024;
    //const height = 768;

    const ctx = canvas.getContext('2d');

    // set background to black

    ctx.fillStyle = 'rgb(0, 0, 0)';

    ctx.fillRect(0, 0, width, height);

    // pick a random colour for each pixel and draw it

    for (let i=0; i<height-1; i++) {

    for (let j=0; j<width-1; j++) {

    // pick a random colour for each point

    // selects random values between 0 and 255 for each r,g,b value

        ctx.fillStyle = `rgb(

            ${Math.floor(Math.random() * 255)},

            ${Math.floor(Math.random() * 255)},

            ${Math.floor(Math.random() * 255)} )`;

    // draw 1x1 point using the selected colour

    // at point (j,i) on the screen

        ctx.fillRect(j, i, 1, 1);

    }

    }
}
