// By: Megan Tibbles
// ray tracing program 

let sphereCenters = [];
let R = [];
let sphereAmbients = [];
let sphereDiffuses = [];
let sphereSpeculars = [];
let sphereShiny = [];

let lightCenter = [];
let lightAmb = [];
let pointLight = [];

let meshTrans = [];
let meshRot = [];
let meshScale = 0;
let meshAmb = [];
let meshDiffuse = [];
let meshSpec = [];
let meshShiny = 0;
let mesh = 0;

let vertices = [];
let textures = []
let normals = [];
let indices = [];

let verticesTemp = [];
let texturesTemp = []
let normalsTemp = [];
let indicesTemp = [];


let numSpheres = 0;

const width = 512;
const height = 512;

function callInit(done){
    if (done == 2){
        if (mesh == 1){
            parse_obj();
        }
        calculate_image();
    }
  }

//reading in needed info from file
function read_file(input){
    let done = 0;
    for (let k = 0; k < input.files.length; k++){
        let file = input.files[k]
        if (file.name.endsWith(".txt")){
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
                        R.push(parseFloat(line[4]));
                        sphereAmbients.push(parseFloat(line[5]), parseFloat(line[6]), parseFloat(line[7]));
                        sphereDiffuses.push(parseFloat(line[8]), parseFloat(line[9]), parseFloat(line[10]));
                        sphereSpeculars.push(parseFloat(line[11]), parseFloat(line[12]), parseFloat(line[13]));
                        sphereShiny.push(parseFloat(line[14]));
                    } else if (type == "mesh"){
                        mesh = 1;
                        meshTrans.push(parseFloat(line[1]), parseFloat(line[2]), parseFloat(line[3]));
                        meshRot.push(parseFloat(line[4]), parseFloat(line[5]), parseFloat(line[6]));
                        meshScale = parseFloat(line[7]);
                        meshAmb.push(parseFloat(line[8]), parseFloat(line[9]), parseFloat(line[10]));
                        meshDiffuse.push(parseFloat(line[11]), parseFloat(line[12]), parseFloat(line[13]));
                        meshSpec.push(parseFloat(line[14]), parseFloat(line[15]), parseFloat(line[16]));
                        meshShiny = parseFloat(line[17]);

                    } else if (type!=''){
                        numSpheres = parseInt(type);
                    }
                    
                }
                console.log(meshTrans, meshRot, meshScale, meshAmb, meshDiffuse, meshSpec, meshShiny);
                done++;
                callInit(done);
                if (input.files.length==1){
                    done++;
                    callInit(done);
                }
            };
        } else if (file.name.endsWith(".obj")){
            let reader = new FileReader();
            reader.readAsText(file);

            reader.onload = function() {
                let type = '';
                let contents = this.result.split('\n');
                //go through each line in file
                for (let i =0; i<contents.length; i++){
                    let line = contents[i].split(' ');
                    type = line[0];
                    line = line.filter(function(con){
                    return con != '';
                    })
                    //go through each value in line
                    for (let j = 1; j < line.length; j++){
                        if (type == "v"){
                            verticesTemp.push(parseFloat(line[j]));
                        } else if (type == "vt"){
                            texturesTemp.push(parseFloat(line[j]));
                        } else if (type == "vn"){
                            normalsTemp.push(parseFloat(line[j]));
                        } else if (type == "f"){
                            indicesTemp.push(line[j]);
                        }
                    }
                }
                done++;
                callInit(done);
            };
        }
    }

}

//find normal
function find_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2 ){
	//let v1 = [x0-x1, y0-y1, z0-z1];
	//let v2 = [x0-x2, y0-y2, z0-z2];
	let v1 = [x1-x0, y1-y0, z1-z0];
	let v2 = [x1-x2, y1-y2, z1-z2];
	let n=[];
	n[0] = v1[1] *v2[2] - v1[2]*v2[1];
	n[1] = v1[2] * v2[0] - v1[0]*v2[2];
	n[2] = v1[0] *v2[1] - v1[1] *v2[0];
	return {x: n[0], y: n[1], z: n[2]};
}

//normalize a vector
function normalize_vector(x, y, z){
    length = Math.sqrt(x**2+y**2+z**2);
    return [x/length, y/length, z/length];
}

//returns the discriminant
function find_discriminant(A, B, C){
    dis = B**2 - (4*A*C);
    return dis;
}

//finds the distance
function find_distance(p1, p2){
    distance = Math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2);
    return distance;
}

//calculates the dot product of two vectors
function dot_product(V1, V2){
    prod = (V1[0]*V2[0]) + (V1[1]*V2[1]) + (V1[2]*V2[2]);
    return prod;
}

//calculates the cross product of two vectors
function cross_product(V1, V2){
    prod= [(V1[2]*V2[1])-(V1[1]*V2[2]), (V1[0]*V2[2])-(V1[2]*V2[0]), (V1[1]*V2[0])-(V1[0]*V2[1])];
    return prod;
}

//calculates the normal of a sphere
function calculate_normal(i, c, Sr, s){
    N = [(i[0] - c[(s*3)])/Sr, (i[1] - c[(s*3)+1])/Sr, (i[2] - c[(s*3)+2])/Sr];
    prd = dot_product(N, [0,0,0]);
    if (prd<0){
        N = [-N[0], -N[1], -N[2]];
    }
    return normalize_vector(N[0], N[1], N[2]);
}

//calculates the light vector
function calculate_light(loc, i){
    L = [(loc[0]-i[0]), (loc[1]-i[1]), (loc[2]-i[2])];
    return normalize_vector(L[0], L[1], L[2]);
}

//calculates the reflection vector
function calculate_reflection(N, L){
    con = 2 * (dot_product(N, L));
    temp = [N[0]*con, N[1]*con, N[2]*con];
    Rvec = [temp[0]-L[0], temp[1]-L[1], temp[2]-L[2]];
    return normalize_vector(Rvec[0], Rvec[1], Rvec[2]);
}

//calculates the viewing vector
function calculate_viewing(i){
    let temp = [0-i[0], 0-i[1], 0-i[2]];
    return normalize_vector(temp[0], temp[1], temp[2]);
}

//calcualtes the colour using the illumination equation
function illumination_equation(Nv, Lv, Rv, V, sAmb, sDiff, sSpec, lAmb, pLight, s, shadow){
    colours = ['r', 'g', 'b'];
    newColour = [];
    let NL = dot_product(Nv, Lv);
    let RV = dot_product(Rv, V);
    if (RV<0){
        RV = 0;
    }
    for (let i = 0; i < colours.length; i++){
        amb = lAmb[i]*sAmb[(s*3)+i];
        let col = 0
        if (shadow){
            col = amb;
        } else {
            diff = pLight[i]*sDiff[(s*3)+i]*NL;
            spec = pLight[i]*sSpec[(s*3)+i]* (RV)**90.0;
            col = amb + diff +spec;
        }
        newColour.push(col*255);
    }
    
    return newColour;
}

//applies the transformations to the object
function transform(points, trans, rot, scale){
    rotPointX = points[0];
    rotPointY = points[1];
    rotPointZ = points[2];
    console.log(rotPointX, rotPointY, rotPointZ);
    for(let x =0; x<points.length-1; x+=3){
        //move to origin
        points[x] -= rotPointX;
        points[x+1] -= rotPointY;
        points[x+2] -= rotPointZ;
        let temp = points[x+1];
        points[x+1] = temp * Math.cos(rot[0]) - points[x+2]*Math.sin(rot[0]);
        points[x+2] = temp * Math.sin(rot[0]) + points[x+2]*Math.cos(rot[0]);
        temp = points[x];
        points[x] = temp * Math.cos(rot[1]) + points[x+2]*Math.sin(rot[1]);
        points[x+2] = temp * Math.sin(rot[1]) + points[x+2]*Math.cos(rot[1]);
        temp = points[x];
        points[x] = temp*Math.cos(rot[2]) - points[x+1]*Math.sin(rot[2]);
        points[x+1] = temp*Math.sin(rot[2]) + points[x+1]*Math.cos(rot[2]);

        points[x] += rotPointX;
        points[x+1] += rotPointY;
        points[x+2] += rotPointZ;
        points[x] = (points[x]*scale)+trans[0];
        points[x+1] = (points[x+1]*scale)+trans[1];
        points[x+2] = (points[x+2]*scale)+trans[2];

    }
    return points;
}

//parses the values from the obj file
function parse_obj(){
    let index = 0;
    let vertexCount = 0;
    verticesTemp = transform(verticesTemp, meshTrans, meshRot, meshScale);
    for (let i = 0; i< indicesTemp.length; i+= 3){
		parts1 = indicesTemp[i].split('/');
		parts2 = indicesTemp[i+1].split('/');
		parts3 = indicesTemp[i+2].split('/');

        vId1 = parseInt(parts1[0]) - 1;
		vId2 = parseInt(parts2[0]) - 1;
		vId3 = parseInt(parts3[0]) - 1;

		tId1 = parseInt(parts1[1]) - 1;
		tId2 = parseInt(parts2[1]) - 1;
		tId3 = parseInt(parts3[1]) - 1; 
        
        vertices.push(verticesTemp[vId1*3], verticesTemp[vId1*3+1], verticesTemp[vId1*3+2]);
		vertices.push(verticesTemp[vId2*3], verticesTemp[vId2*3+1], verticesTemp[vId2*3+2]);
		vertices.push(verticesTemp[vId3*3], verticesTemp[vId3*3+1], verticesTemp[vId3*3+2]);

		textures.push(texturesTemp[tId1*2], 1-texturesTemp[tId1*2+1]);
		textures.push(texturesTemp[tId2*2], 1-texturesTemp[tId2*2+1]);
		textures.push(texturesTemp[tId3*2], 1-texturesTemp[tId3*2+1]);
		//if normals are included
		if (parts1.length > 2){
			nId1 = parseInt(parts1[2]) - 1;
			nId2 = parseInt(parts2[2]) - 1;
			nId3 = parseInt(parts3[2]) - 1;
            normals.push(normalsTemp[nId1*3], normalsTemp[nId1*3+1], normalsTemp[nId1*3+2]);
			normals.push(normalsTemp[nId2*3], normalsTemp[nId2*3+1], normalsTemp[nId2*3+2]);
			normals.push(normalsTemp[nId3*3], normalsTemp[nId3*3+1], normalsTemp[nId3*3+2]);


		} else {
			// if normals must be calculated
			let n = find_normal(verticesTemp[vId1*3], verticesTemp[vId1*3+1], verticesTemp[vId1*3+2], verticesTemp[vId2*3], verticesTemp[vId2*3+1], verticesTemp[vId2*3+2], verticesTemp[vId3*3], verticesTemp[vId3*3+1], verticesTemp[vId3*3+2]);
			
            length = Math.sqrt(n.x**2 + n.y**2 + n.z**2);
			normals.push(n.x/length, n.y/length, n.z/length);
			normals.push(n.x/length, n.y/length, n.z/length);
			normals.push(n.x/length, n.y/length, n.z/length);
		}
		indices.push(index, index+1, index+2);
		index+=3;

	}
	vertexCount = vertices.length/3;
    console.log(vertices, indices, textures, normals);
}

//check if a ray intersections with an object. Returns intersection point and the normal vector
function check_intersection(i, V, vp){
    epsilon = 0.0000001;
    i1 = indices[i];
    i2 = indices[i+1];
    i3 = indices[i+2];
    v1=[vertices[(i1*3)],vertices[(i1*3)+1], vertices[(i1*3)+2]];
    v2=[vertices[(i2*3)],vertices[(i2*3)+1], vertices[(i2*3)+2]];
    v3=[vertices[(i3*3)],vertices[(i3*3)+1], vertices[(i3*3)+2]];

    edge1 = [v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]];
    edge2 = [v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2]];
    Ve2 = cross_product(V, edge2);
    let det = dot_product(edge1, Ve2);
    
    if (det>-epsilon && det <epsilon){
        return [0];
    }
    invDet = 1.0/det;
    let S = [vp[0]-v1[0], vp[1]-v1[1], vp[2]-v1[2]];
    u = invDet*dot_product(S, Ve2);

    if (u < 0.0 || u > 1.0){
        return [0];
    }

    Se1 = cross_product(S, edge1);
    v = invDet*dot_product(V, Se1);
    
    if (v < 0 || (u+v)>1){
        return [0];
    }
    t = invDet*dot_product(edge2, Se1);
                        
    if (t>epsilon){
        Nv = [-normals[i1*3], normals[(i1*3)+1], normals[(i1*3)+2]];
        temp = [V[0]*t, V[1]*t, V[2]*t];
        ri = [vp[0]+temp[0], vp[1]+temp[1], vp[1]+temp[2]];
        return [ri, Nv];
    }
    return [0];
}

//check if the ray intersects with the sphere
function check_sphere_intersection(s, V, vp){
    let ri = []
    let A = 1.0
    let B = 2 * ( (V[0] * (vp[0] - sphereCenters[(s*3)])) + (V[1] * (vp[1] - sphereCenters[(s*3)+1])) + (V[2] * (vp[2] - sphereCenters[(s*3)+2])));
    let C  = (vp[0] - sphereCenters[(s*3)])**2 + (vp[1] - sphereCenters[(s*3)+1])**2 + (vp[2] - sphereCenters[(s*3)+2])**2 - R[s]**2;
    dis = find_discriminant(A, B, C);

    //only intersects once
    if (dis < 0.01 && dis > -0.01){
        t0 = (-B - Math.sqrt(Math.abs(dis)))/(2*A);
        if (t0<0){
            return [];
        }
        ri = [(vp[0] + V[0]*t0), (vp[1] + V[1]*t0), (vp[2] + V[2]*t0)];

    } else if (dis > 0) {   //intersects twice
        t0 = (-B - Math.sqrt(dis))/(2*A);
        t1 = (-B + Math.sqrt(dis))/(2*A);
        if (t0>=0 && t1>=0){
            let ri0 = [(vp[0] + V[0]*t0), (vp[1] + V[1]*t0), (vp[2] + V[2]*t0)];
            let ri1 = [(vp[0] + V[0]*t1), (vp[1] + V[1]*t1), (vp[2] + V[2]*t1)];
            //save closest intersection value
            distance1 = find_distance(vp, ri0);
            distance2 = find_distance(vp, ri1);
            if (distance1<distance2){
                ri = ri0;
            } else {
                ri = ri1;
            }
        //if < 0, behind the origin point
        } else if (t0>=0){
            ri = [(vp[0] + V[0]*t0), (vp[1] + V[1]*t0), (vp[2] + V[2]*t0)];
        } else if (t1>=0){
            ri = [(vp[0] + V[0]*t1), (vp[1] + V[1]*t1), (vp[2] + V[2]*t1)];
        }
    } 
    return ri;
}



//viewing plane where z = -1
function calculate_image(){
    const canvas = document.querySelector('.myCanvas');
    const ctx = canvas.getContext('2d');

    ctx.fillStyle = 'rgb(0, 0, 0)';
    ctx.fillRect(0, 0, width, height);

    let z = -1
    let vp = [0.0,0.0,-1.0];
    let w2 = 2/width;
    let h2 = 2/height;
    let x, y;
    
    for (x=0; x<width-1; x+=1){
        for (y=0; y<height-1;y+=1){
            //calculate a ray from the viewpoint (0,0,0) for this pixel.
            //Normalize it
            let V = [((x*w2)-1), 0-((y*h2)-1), z];
            V = normalize_vector(V[0], V[1], V[2]);
            let closestDist = 100;
            //calculate if the ray intersects the sphere
            for (s=0; s<numSpheres; s++){
                let ri = check_sphere_intersection(s, V, vp)
                dist = find_distance(vp, ri);  
                let shadow = false;
                let blendColour = [];
                if (ri.length == 3){
                    Vv = calculate_viewing(ri);
                    Nv = calculate_normal(ri, sphereCenters, R[s], s);
                    Lv = calculate_light(lightCenter, ri);
                    Rv = calculate_reflection(Nv, Lv);
                    for (s2=0; s2<numSpheres; s2++){
                        if (s2!=s){
                            let ri2 = check_sphere_intersection(s2, Lv, ri);
                            if (ri2.length==3){
                                shadow = true;
                            }

                            let RvRi = check_sphere_intersection(s2, Rv, ri);
                            if (RvRi.length==3){
                                RvVv = calculate_viewing(RvRi);
                                RvNv = calculate_normal(RvRi, sphereCenters, R[s2], s2);
                                RvLv = calculate_light(lightCenter, RvRi);
                                RvRv = calculate_reflection(RvNv, RvLv); 
                                blendColour = illumination_equation(Nv, Lv, Rv, Vv, sphereAmbients, sphereDiffuses, sphereSpeculars, lightAmb, pointLight, s2, false);
                            }
                        }
                    }
                    for (let i2=0; i2<indices.length; i2+=3){
                        let ri3 = check_intersection(i2, Lv, ri);
                        if (ri3[0]!=0){
                            shadow = true;
                        }
                        let RvRi = check_intersection(i2, Rv, ri);
                        if (RvRi.length!=1){
                            RvVv = calculate_viewing(RvRi[0]);
                            RvNv = RvRi[1]
                            RvLv = calculate_light(lightCenter, RvRi[0]);
                            RvRv = calculate_reflection(RvNv, RvLv); 
                            blendColour = illumination_equation(RvNv, RvLv, RvRv, RvVv, meshAmb, meshDiffuse, meshSpec, lightAmb, pointLight, 0, false);
                        }
                    }
                    if (dist<closestDist){
                        colour = illumination_equation(Nv, Lv, Rv, Vv, sphereAmbients, sphereDiffuses, sphereSpeculars, lightAmb, pointLight, s, shadow);
                        if (blendColour.length !=0){
                            colour = [(1.0 -sphereShiny[s])*colour[0]+ (sphereShiny[s])*blendColour[0], (1.0 -sphereShiny[s])*colour[1] + (sphereShiny[s])*blendColour[1], (1.0 -sphereShiny[s])*colour[2] + (sphereShiny[s])*blendColour[2]];
                        }
                        rgb = `rgb(`+colour[0]+`,`+colour[1]+`,`+colour[2]+`)`;
                        ctx.fillStyle =  rgb;
                        ctx.fillRect(x, y, 1, 1);
                    } 
                    closestDist = dist;
                } 
            }
            if (mesh == 1){
                let ri = [];
                for (let i=0; i<indices.length; i+=3){
                    let rgb =0;
                    let shadow = false;
                    let blendColour = [];
                    ri = check_intersection(i, V, vp);
                    if (ri.length!=1){
                        Nv = ri[1]
                        Vv = calculate_viewing(ri[0]);
                        Lv = calculate_light(lightCenter, ri[0]);
                        Rv = calculate_reflection(Nv, Lv);
                        for (let i2=0; i2<indices.length; i2+=3){
                            if (i2!=i){
                                let ri2 = check_intersection(i2, Lv, ri[0]);
                                if (ri2[0]!=0){
                                    shadow = true;
                                }
                                let RvRi = check_intersection(i2, Rv, ri[0]);
                                if (RvRi.length!=1){
                                    RvVv = calculate_viewing(RvRi[0]);
                                    RvNv = RvRi[1]
                                    RvLv = calculate_light(lightCenter, RvRi[0]);
                                    RvRv = calculate_reflection(RvNv, RvLv); 
                                    blendColour = illumination_equation(RvNv, RvLv, RvRv, RvVv, meshAmb, meshDiffuse, meshSpec, lightAmb, pointLight, 0, false);
                                }
                            }
                        }
                        for (s2=0; s2<numSpheres; s2++){
                            let ri2 = check_sphere_intersection(s2, Lv, ri[0]);
                            if (ri2.length==3){
                                shadow = true;
                                
                            }
                            let reflectRay = [-Rv[0], Rv[1],Rv[2]];
                            let RvRi = check_sphere_intersection(s2, reflectRay, ri[0]);
                            if (RvRi.length==3){
                                //console.log(x,y);
                                //RvVv = calculate_viewing(RvRi);
                                //RvNv = calculate_normal(RvRi, sphereCenters, R[s2], s2);
                                //RvLv = calculate_light(lightCenter, RvRi);
                                //RvRv = calculate_reflection(RvNv, RvLv); 
                                blendColour = illumination_equation(Nv, Lv, Rv, Vv, sphereAmbients, sphereDiffuses, sphereSpeculars, lightAmb, pointLight, s2, false);
                            }
                        }
                        let colour = illumination_equation(Nv, Lv, Rv, Vv, meshAmb, meshDiffuse, meshSpec, lightAmb, pointLight, 0, shadow);
                        if (blendColour.length !=0){
                            colour = [(1.0 -meshShiny)*colour[0]+ (meshShiny)*blendColour[0], (1.0 -meshShiny)*colour[1] + (meshShiny)*blendColour[1], (1.0 -meshShiny)*colour[2] + (meshShiny)*blendColour[2]];
                        }
                        rgb1 = `rgb(`+colour[0]+`,`+colour[1]+`,`+colour[2]+`)`;
                        dist = find_distance(vp, ri[0]);
                        if (dist<closestDist){
                            rgb = rgb1;
                            closestDist = dist;
                        }
                    }
                    if (rgb!=0){
                        ctx.fillStyle =  rgb;
                        ctx.fillRect(x, y, 1, 1);
                    }
                }
            }

        }
    }

}
