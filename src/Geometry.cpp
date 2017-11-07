#include <fstream>
#include <iostream>
#include <algorithm>

#include "TPolyLine.h"
#include "TVector2.h"


#ifdef STANDALONE
#include "Geometry.h"
#include "json/json.h"
#else
#include "HGCalSimulation/FastShower/interface/Geometry.h"
#include "HGCalSimulation/FastShower/interface/json/json.h"
#endif


using namespace std;


Geometry::~Geometry(){
}


string fixedLength(int value, int digits = 5, string type="FH") {
    //  In the ex. below wants
    //  for int value 0 -> it shoud display 000
    //  for value -9 -> it should return -009
    //  for value say -50 -> it should return -050
    //  for value say -110 -> it should return -110

    if (type != "FH" && type != "VHH") {
        cout << "Unknown cell type name " << type << endl;
        return "";
    }

    unsigned int uvalue = value;
    string result;
    while (digits-- > 0) {
        result += ('0' + uvalue % 10);
        uvalue /= 10;
    }
    if (value < 0) {
        result += '-';
    }
    reverse(result.begin(), result.end());
    result = type + result;
    return result;
}


void Geometry::setLayer(int klayer) {
    if (klayer<0 || unsigned(klayer)>=parameters_.layers_z.size()) {
        stringstream error;
        error << "[Geometry] error, invalid klayer " << klayer;
        throw error.str();
    }
    klayer_ = klayer;
}


std::array<double, 4> Geometry::dimensions(double side) {
    std::array<double, 4> dim;

    dim[0] = side*sqrt(3.); //asqrt3_
    dim[1] = dim[0]/2.;  //asqrt3over2_
    dim[2] = side/2.;          //aover2_
    dim[3] = dim[2]*3.;      //a3over2

    return dim;
}

std::array< std::array<double, 6>, 4> Geometry::hexagonoffset(double side) {

    // pic on top
    std::array<double, 4> dim = dimensions(side);

    return {{ {{dim[1], dim[1], 0., -dim[1],-dim[1], 0}},
              {{-dim[2], dim[2], side, dim[2], -dim[2],-side}}
           }};

    // Hexagons flat side on top
    // const array<double, nvertices> hexagonoffsetx = {{aover2_,a_,aover2_,-aover2_,-a,-aover2_}};
    // const array<double, nvertices> hexagonoffsety = {{-asqrt3over2_,0,asqrt3over2_,asqrt3over2_,0.,-asqrt3over2_}};

    // return offset;
}

std::array< std::array<double, 6>, 4> Geometry::triangleoffset(double side){

    // up, up, down, down
    std::array<double, 4> dim = dimensions(side);
    return {{ {{dim[2], 0., -dim[2]}},
              {{-dim[1]/3., dim[0]/3., -dim[1]/3.}},
              {{dim[2], -dim[2], 0.}},
              {{dim[1]/3., dim[1]/3., -dim[0]/3.}},
           }};
}

std::array<double, 3> Geometry::derivative(double side, Parameters::Geometry::Type itype){
    std::array<double, 3> derivative;

    std::array<double, 4> dim = dimensions(side);
    // partial derivatives of x,y vs i,j
    switch(itype) {
        case Parameters::Geometry::Type::Hexagons: {
            // peak on top
            derivative[0] = dim[0];
            derivative[1] = dim[1];
            derivative[2] = dim[3];

            // // flat side on top 
            // dydj = asqrt3_;
            // dxdj = asqrt3over2_;
            // dxdi = a3over2_;
            break;
        }
        case Parameters::Geometry::Type::Triangles: {
            derivative[0] = dim[2];
            derivative[1] = dim[2];
            derivative[2] = dim[1];
            break;
        }
        default:
            break;
    };
    return derivative;
}

std::array< std::array<double, 5>, 2> Geometry::dxdyFirstZone(std::array<double, 9> xs, std::array<double, 9> ys) {
    std::array<double, 5> dx;
    std::array<double, 5> dy;

    for(unsigned i=0; i<5; i++)
        dx[i] = xs[i+1]-xs[0];
    for(unsigned i=0; i<5; i++)
        dy[i] = ys[i+1]-ys[0];

    return {{ dx, dy}};
}

std::array< std::array<double, 5>, 2> Geometry::dxdySecondZone(std::array<double, 9> xs, std::array<double, 9> ys) {
    std::array<double, 5> dx;
    std::array<double, 5> dy;

    dx[0] = xs[2]-xs[3];
    dx[1] = xs[6]-xs[3];
    dx[2] = xs[7]-xs[3];
    dx[3] = xs[5]-xs[3];
    dx[4] = xs[8]-xs[3];

    dy[0] = ys[2]-ys[3];
    dy[1] = ys[6]-ys[3];
    dy[2] = ys[7]-ys[3];
    dy[3] = ys[5]-ys[3];
    dy[4] = ys[8]-ys[3];

    return {{ dx, dy}};
}

std::array< int, 4> Geometry::ijWindows(int layer_region, std::array<double, 9> xs, std::array<double, 9> ys, double side, Parameters::Geometry::Type itype) {
    // compute i,j window needed to cover the x,y window

    std::array< int, 4> windows;    //imin, imax, jmin, jmax

    std::array<double, 3> der = derivative(side, itype);
    std::array< std::array<double, 5>, 2> dxdy;
    array<double,6> js_first;
    array<double,6> is_first;
    array<double,6> js_second;
    array<double,6> is_second;

    // layer_region = part of layer : 1 = small cells part, 2 = large cells part
    if (layer_region == 1) {
        dxdy = dxdyFirstZone(xs, ys);
        js_first[0] = 0;
        is_first[0] = 0;
        for(unsigned i=1; i<js_first.size(); i++) js_first[i] = dxdy[1][i-1]/der[2];
        for(unsigned i=1; i<is_first.size(); i++) is_first[i] = (dxdy[0][i-1]-js_first[i]*der[1])/der[0];

        windows[0] = round(*(min_element(is_first.begin(), is_first.end())));
        windows[2] = round(*(min_element(js_first.begin(), js_first.end())));
        windows[1] = round(*(max_element(is_first.begin(), is_first.end())));
        windows[3] = round(*(max_element(js_first.begin(), js_first.end())));
    }
    else {
        dxdy = dxdySecondZone(xs, ys);
        js_second[0] = 0;
        is_second[0] = 0;
        for(unsigned i=1; i<js_second.size(); i++) js_second[i] = dxdy[1][i-1]/der[2];
        for(unsigned i=1; i<is_second.size(); i++) is_second[i] = (dxdy[0][i-1]-js_second[i]*der[1])/der[0];

        windows[0] = round(*(min_element(is_second.begin(), is_second.end())));
        windows[2] = round(*(min_element(js_second.begin(), js_second.end())));
        windows[1] = round(*(max_element(is_second.begin(), is_second.end())));
        windows[3] = round(*(max_element(js_second.begin(), js_second.end())));
    }
    return windows;
}


std::array< double, 4> Geometry::XYrPhi(int i, int j, double side, Parameters::Geometry::Type itype, std::array<double, 9> xs, std::array<double, 9> ys, double layer_region) {
    std::array< double, 4> xyrPhi;

    std::array<double, 3> der = derivative(side, itype);
    // peak on top
    if (layer_region == 1) {
        xyrPhi[0] = xs[0] + i*der[0] + j*der[1];
        xyrPhi[1] = ys[0] + j*der[2];
    }
    else {
        xyrPhi[0] = xs[3] + i*der[0] + j*der[1];
        xyrPhi[1] = ys[3] + j*der[2];
    }

    // //flat side on top
    // double x = xs[0] + i*dxdi;
    // double y = ys[0] + j*dydj + i*dxdj;

    // up and down triangle barycenters are not aligned
    if(itype==Parameters::Geometry::Type::Triangles && i%2) {
        std::array<double, 4> dim = dimensions(side);
        xyrPhi[1] += dim[0]/6.;
    }

    xyrPhi[2] = sqrt(xyrPhi[0]*xyrPhi[0] + xyrPhi[1]*xyrPhi[1]);
    xyrPhi[3] = copysign(acos(xyrPhi[0]/xyrPhi[2]), xyrPhi[1]);

    return xyrPhi;
}


void Geometry::constructFromJson(bool debug, int layer_id) {

    itype_ = Parameters::Geometry::Type::External;
    const string& filename = parameters_.file;

    // first set klayer and zlayer according to user parameters
    setLayer(layer_id);
    double zlayer = parameters_.layers_z[klayer_]; // else offset from the layer z position
    setZlayer(zlayer);

    // json format from Marina 06/2016
    // units are mm, converted in cm
    Json::Reader reader;
    Json::Value obj;
    ifstream ifs(filename);
    //bool success = reader.parse(ifs, obj);
    reader.parse(ifs, obj);
    const Json::Value& version = obj["Version"]; 
    if(debug) {
        cout << " " << endl;
        cout << "Reading geometry from JSON file: " << filename << endl;
        cout << "Geometry version " << version.asString() << endl;
    }

    // module meta data
    const Json::Value& module = obj["Module"];
    if(module.isNull())
        throw string("No module information found in json file");
    const Json::Value& module_area = module["Module_area"]; 
    const Json::Value& module_cells_and_half_cells = module["Module_cells_(1/1, 1/2)"]; 
    const Json::Value& module_full_cells = module["Module_cells_1/1"]; 
    const Json::Value& module_half_cells = module["Module_cells_1/2"]; 
    const Json::Value& module_third_cells = module["Module_cells_1/3"]; 
    const Json::Value& module_center_coord = module["Module_center_coord"]; 
    const Json::Value& module_orientation = module["Module_orientation"]; 
    const Json::Value& module_vertices = module["Module_vertex_coord"];
    //bool module_vertices_ok = false;
    vector<pair<double,double>> module_vertex_coordinates;
    if(!module_vertices.isNull() && module_vertices.isArray()) {
        //module_vertices_ok = true;
        for (unsigned i=0; i<module_vertices.size(); i++) {
            const Json::Value& coord = module_vertices[i]; 
            if(coord.isNull() || !coord.isArray() || coord.size()!=2) {
                //module_vertices_ok = false;
                break;
            }
            module_vertex_coordinates.emplace_back(coord[0].asDouble(), coord[1].asDouble());
        }
    }
    if(!module_cells_and_half_cells.isNull() && !module_full_cells.isNull() && !module_half_cells.isNull()) {
        if(module_cells_and_half_cells.asUInt()!=module_half_cells.asUInt()+module_full_cells.asUInt()) {
            cout<<"Inconsistency in the number of full and half cells\n";
        }
    }
    // FIXME: should add more checks
    if(debug) {
        cout << "Module area : " << module_area.asDouble() << endl;
        cout << "Module cells (1/1, 1/2) : " << module_cells_and_half_cells.asUInt() << endl;
        cout << "Module cells_1/1 : " << module_full_cells.asUInt() << endl;
        cout << "Module cells 1/2 : " << module_half_cells.asUInt() << endl;
        cout << "Module cells 1/3 : " << module_third_cells.asUInt() << endl;
        cout << "Module center coord : " << module_center_coord[0].asDouble() << 
          " " << module_center_coord[1].asDouble() << endl;
        cout << "Module orientation : " << module_orientation.asDouble() << endl;
        cout << "Module vertex coord : " << endl;
        cout << " " << endl;
    }

    // full hexagon cells meta data
    const Json::Value& hexagons = module["FH"];
    if(hexagons.isNull() || !hexagons.isObject())
        throw string("Cannot find full hexagons information");
    const Json::Value& full_hexagons_area = hexagons["FH_area"]; 
    const Json::Value& full_hexagons_area_module = hexagons["FH_area_module"]; 
    const Json::Value& full_hexagons_count = hexagons["FH_count"]; 
    // FIXME: should add more checks
    if(debug) {
        cout << "Full hexagon area " << full_hexagons_area.asDouble() << endl;
        cout << "Full hexagon area module " << full_hexagons_area_module.asDouble() << endl;
        cout << "Full hexagons nbr of cells " << full_hexagons_count.asUInt() << endl;
        cout << " " << endl;
    }

    // construct full hexagon cells
    for (unsigned int icell=0; icell<full_hexagons_count.asUInt(); icell++) {
        string cell_name = fixedLength(icell,5,"FH");
        // FIXME: add more checks
        const Json::Value& hexagon_attributes = hexagons[cell_name]; 
        int i_index = hexagon_attributes[0]["mapping_coord"][0].asInt();
        int j_index = hexagon_attributes[0]["mapping_coord"][1].asInt();

        TVectorD position(3);
        position(0) = hexagon_attributes[1]["center_coord"][0].asDouble()/10.; // cm
        position(1) = hexagon_attributes[1]["center_coord"][1].asDouble()/10.; // cm
        position(2) = getZlayer();

        double orientation = hexagon_attributes[2]["orientation"].asDouble();
        vector<TVectorD> vertices;
        for (unsigned i=0; i<hexagon_attributes[3]["vertex_coord_abs"].size(); i++) {
            vertices.emplace_back(2);
            vertices.back()(0) = hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble()/10.; // cm
            vertices.back()(1) = hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble()/10.; // cm
        }
        cells_.emplace(Cell::id(i_index, j_index,klayer_), Cell(move(position), move(vertices), orientation, i_index, j_index,klayer_));

        if(debug) {
            cout << "New cell " << cell_name << " : " << endl;
            cout << " mapping coordinates : " << 
              hexagon_attributes[0]["mapping_coord"][0].asInt() << " " <<
              hexagon_attributes[0]["mapping_coord"][1].asInt() << endl;
            cout << " center coordinates : " << 
              hexagon_attributes[1]["center_coord"][0].asDouble() << " " <<
              hexagon_attributes[1]["center_coord"][1].asDouble() << endl;
            cout << " orientation : " << 
              hexagon_attributes[2]["orientation"].asDouble() << endl;
            cout << " vertex_coordinates : " << endl;
            for(const auto& vertex : vertices) {
                cout << "  vertex " << vertex(0) << " " << vertex(1) << endl;
            }
        }
    }
    if(debug) cout << " " << endl;

    // half hexagon cells meta data
    const Json::Value& half_hexagons = module["Edge_VHH"];
    const Json::Value& half_hexagons_area = half_hexagons["VHH_area"];
    const Json::Value& half_hexagons_area_module = half_hexagons["VHH_area_module"];
    const Json::Value& half_hexagons_count = half_hexagons["VHH_count"];
    // FIXME: should add more checks
    if(debug) {
        cout << "Half hexagon area (edges) : " << half_hexagons_area.asDouble() << endl;
        cout << "Half hexagon area module (edges) : " << half_hexagons_area_module.asDouble() << endl;
        cout << "Half hexagons nbr of cells (edges) : " << half_hexagons_count.asUInt() << endl;
        cout << " " << endl;
    }

    // construct half hexagon cells
    for (unsigned int icell=0; icell<half_hexagons_count.asUInt(); icell++) {
        string cell_name = fixedLength(icell,5,"VHH");
        // FIXME: add more checks
        const Json::Value& hexagon_attributes = half_hexagons[cell_name];
        int i_index = hexagon_attributes[0]["mapping_coord"][0].asInt();
        int j_index = hexagon_attributes[0]["mapping_coord"][1].asInt();

        TVectorD position(3);
        position(0) = hexagon_attributes[1]["center_coord"][0].asDouble()/10.; // cm
        position(1) = hexagon_attributes[1]["center_coord"][1].asDouble()/10.; // cm
        position(2) = getZlayer();

        double orientation = hexagon_attributes[2]["orientation"].asDouble();
        vector<TVectorD> vertices;

        for (unsigned i=0; i<hexagon_attributes[3]["vertex_coord_abs"].size(); i++) {
            vertices.emplace_back(2);
            vertices.back()(0) = hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble()/10.;// cm
            vertices.back()(1) = hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble()/10.;// cm
        }
        cells_.emplace(Cell::id(i_index, j_index,klayer_),Cell(move(position), move(vertices), orientation, i_index, j_index,klayer_));

        if(debug) {
            cout << "New cell " << cell_name << " : " << endl;
            cout << " mapping coordinates : " <<
              hexagon_attributes[0]["mapping_coord"][0].asInt() << " " <<
              hexagon_attributes[0]["mapping_coord"][1].asInt() << endl;
            cout << " center coordinates : " <<
              hexagon_attributes[1]["center_coord"][0].asDouble() << " " <<
              hexagon_attributes[1]["center_coord"][1].asDouble() << endl;

            cout << " orientation : " << hexagon_attributes[2]["orientation"].asDouble() << endl;
            cout << " vertex_coordinates : " << endl;
            for(const auto& vertex : vertices) {
              cout << "  vertex " << vertex(0) << " " << vertex(1) << endl;
            }
        }
    }
    if(debug) cout << " " << endl;
}


void Geometry::constructFromParameters(bool debug, int layer_id, int display_layer) {
  // a tesselation of the plane with polygons

    Parameters::Geometry::Type itype(parameters_.type);
    double zlayer = parameters_.layers_z[layer_id]; // else offset from the layer z position
    setZlayer(zlayer);
    int klayer(layer_id);
    setLayer(klayer);

    if(debug) {
        cout << " " << endl;
        cout << "Building parametrized geometry :\n";
        cout << "  " << parameters_.eta_min << " < eta < " << parameters_.eta_max << "\n";
        cout << "  " << parameters_.phi_min << " < phi < " << parameters_.phi_max << "\n";
        if (itype==Parameters::Geometry::Type::Hexagons)
            cout << "with hexagonal cells " << endl;
        else if (itype==Parameters::Geometry::Type::Triangles)
            cout << "with triangular cells " << endl;
    }

    // compute x,y positions of the global geometry window (in the eta-phi region)
    double theta_min = 2.*atan(exp(-parameters_.eta_max));
    double theta_max = 2.*atan(exp(-parameters_.eta_min));
    double r_min = zlayer*tan(theta_min);
    double r_max = zlayer*tan(theta_max);
    double phi_min = parameters_.phi_min;
    double phi_max = parameters_.phi_max;

    // Define sub-zone in in this geometry windows : 100um, 200um, 300um
    double limit_first_zone = parameters_.limit_first_zone;

    double x1 = r_min*cos(phi_max);
    double y1 = r_min*sin(phi_max);
    double x4 = r_min*cos((phi_min+phi_max)/2.);
    double y4 = r_min*sin((phi_min+phi_max)/2.);

    // Compute x,y positions of the sub-zone in layer (100,200,300um)
    double denom = (2*x4*x4-((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)));
    double xs_sup_first = (limit_first_zone/(2*x4*x4))*denom;
    double ys_sup_first = sqrt((limit_first_zone*limit_first_zone)-(xs_sup_first *xs_sup_first));

    std::array<double, 9> xs = {{
        r_min*cos(phi_min),
        x1,
        xs_sup_first,
        xs_sup_first,
        x4,
        limit_first_zone,
        r_max*cos(phi_max),
        r_max*cos(phi_min),
        r_max*cos((phi_min+phi_max)/2.)
    }};

    std::array<double, 9> ys = {{
        r_min*sin(phi_min),
        y1,
        ys_sup_first,
        -ys_sup_first,
        y4,
        y4,
        r_max*sin(phi_max),
        r_max*sin(phi_min),
        r_max*sin((phi_min+phi_max)/2.)
    }};

    //For the output geometry TH1
    double x_min = numeric_limits<double>::max();
    double x_max = numeric_limits<double>::lowest();
    double y_min = numeric_limits<double>::max();
    double y_max = numeric_limits<double>::lowest();

    int n_vertices;
    if (itype==Parameters::Geometry::Type::Hexagons)
        n_vertices = 6;
    else
        n_vertices = 3;

    double side1 = parameters_.small_cell_side;

    //cell_size_region = 1 -> small cell part
    //cell_size_region = 2 -> large cell part
    int cell_size_region = 1;

    std::array< std::array<double, 6>, 4> offset1;
    if (itype==Parameters::Geometry::Type::Hexagons) {
        offset1 = hexagonoffset(side1);
    }
    else
        offset1 = triangleoffset(side1);

    std::array< int, 4> windows1 = ijWindows(cell_size_region, xs, ys, side1, itype);


    double orientation = 90.;

    // build cells inside the requested window
    int real_i = 0;
    int real_j = 0;
    for (int i = windows1[0]; i <= windows1[1]; i++) {
        real_j = 0;
        for (int j = windows1[2]; j <= windows1[3]; j++) {

            std::array< double, 4> xyrPhi = XYrPhi(i, j, side1, itype, xs, ys,cell_size_region);

            // check if cell is inside boundaries. If not, skip it
            if(!(xyrPhi[2] >= r_min &&
                 xyrPhi[2] <= limit_first_zone &&
                 TVector2::Phi_mpi_pi(xyrPhi[3]-phi_min)>=0 &&
                 TVector2::Phi_mpi_pi(xyrPhi[3]-phi_max)<=0)) {
                continue;
            }

            if(xyrPhi[0]<x_min) x_min = xyrPhi[0];

            if (debug) {
                cout << "Creating new cell of type "
                     << static_cast<underlying_type<Parameters::Geometry::Type>::type>(itype)
                     << " : " << endl;
                cout << " mapping coordinates : " << i << " " << j << endl;
            }

            TVectorD position(3);
            position(0) = xyrPhi[0];
            position(1) = xyrPhi[1];
            position(2) = zlayer;

            if (debug) {
                cout << " center coordinates : " << xyrPhi[0] << " " << xyrPhi[1] << " "<< zlayer << endl;
            }

            if (itype == Parameters::Geometry::Type::Triangles && i%2 != 0) orientation =  -90.; // for downward triangles

            vector<TVectorD> vertices;
            for (int iv=0; iv<n_vertices; iv++) {
                vertices.emplace_back(2);
                switch(itype) {
                  case Parameters::Geometry::Type::Hexagons: {
                      vertices.back()(0) = xyrPhi[0] + offset1[0][iv];
                      vertices.back()(1) = xyrPhi[1] + offset1[1][iv];
                      break;
                  }
                  case Parameters::Geometry::Type::Triangles:
                  {
                    switch(i%2) {
                      case 0:
                      {
                        vertices.back()(0) = xyrPhi[0] + offset1[0][iv];
                        vertices.back()(1) = xyrPhi[1] + offset1[1][iv];
                        break;
                      }
                      default:
                      {
                        vertices.back()(0) = xyrPhi[0] + offset1[2][iv];
                        vertices.back()(1) = xyrPhi[1] + offset1[3][iv];
                        break;
                      }
                    }
                    break;
                  }
                  default:
                      break;
                };
                if (debug) {
                  cout<<"vertex "<<iv<<" "<<vertices.back()(0)<<" "<<vertices.back()(1)<<endl;
                }
            }

            auto c = cells_.emplace(
                Cell::id(real_i, real_j, layer_id),
                Cell(std::move(position), std::move(vertices), orientation, real_i, real_j, layer_id));
            if(!c.second) {
                cout << "Warning: Cell with indices" << real_i << " "
                        << real_j << ", "<< layer_id
                        << " already exists (id="<<Cell::id(real_i, real_j, layer_id) << ")\n";
                cout << "Warning: This may indicate a bug in the id definition\n";
            }
            real_j++;
        }
        real_i++;
    }

    double side2 = parameters_.large_cell_side;
    cell_size_region = 2;

    std::array< std::array<double, 6>, 4> offset2;
    if (itype==Parameters::Geometry::Type::Hexagons)
        offset2 = hexagonoffset(side2);
    else
        offset2 = triangleoffset(side2);

    std::array< int, 4> windows2 = ijWindows(cell_size_region, xs, ys, side2, itype);

  // build cells inside the requested window
    int real_i2 = 0;
    int real_j2 = 0;
    for (int i=windows2[0]; i<=windows2[1]; i++) {
        real_j2 = 0;
        for (int j = windows2[2]; j<=windows2[3]; j++) {

            std::array< double, 4> xyrPhi = XYrPhi(i, j, side2, itype, xs, ys, cell_size_region);

            // check if cell is inside boundaries. If not, skip it
            if(!(xyrPhi[2]>=limit_first_zone &&
                 xyrPhi[2]<=r_max &&
                 TVector2::Phi_mpi_pi(xyrPhi[3]-phi_min)>=0 &&
                 TVector2::Phi_mpi_pi(xyrPhi[3]-phi_max)<=0)) {
                continue;
            }

            if(xyrPhi[0]>x_max) x_max = xyrPhi[0];
            if(xyrPhi[1]>y_max) y_max = xyrPhi[1];
            if(xyrPhi[1]<y_min) y_min = xyrPhi[1];

            if (debug) {
              cout << "Creating new cell of type "
                   << static_cast<underlying_type<Parameters::Geometry::Type>::type>(itype)
                   << " : " << endl;
              cout << " mapping coordinates : " << i << " " << j << endl;
            }

            TVectorD position(3);
            position(0) = xyrPhi[0];
            position(1) = xyrPhi[1];
            position(2) = zlayer;
            if (debug) {
                cout << " center coordinates : " << xyrPhi[0] << " " << xyrPhi[1] << endl;
            }

            if (itype == Parameters::Geometry::Type::Triangles && i%2 != 0) orientation =  -90.; // for downward triangles

            vector<TVectorD> vertices;
            for (int iv=0; iv<n_vertices; iv++) {
                vertices.emplace_back(2);
                switch(itype) {
                  case Parameters::Geometry::Type::Hexagons: {
                      vertices.back()(0) = xyrPhi[0]+offset2[0][iv];
                      vertices.back()(1) = xyrPhi[1]+ offset2[1][iv];
                      break;
                  }
                  case Parameters::Geometry::Type::Triangles:
                  {
                    switch(i%2) {
                      case 0:
                      {
                        vertices.back()(0) = xyrPhi[0] + offset2[0][iv];
                        vertices.back()(1) = xyrPhi[1] + offset2[1][iv];
                        break;
                      }
                      default:
                      {
                        vertices.back()(0) = xyrPhi[0] + offset2[2][iv];
                        vertices.back()(1) = xyrPhi[1] + offset2[3][iv];
                        break;
                      }
                    }
                    break;
                  }
                  default:
                      break;
                };
                if (debug) {
                    cout<<"vertex "<<iv<<" "<<vertices.back()(0)<<" "<<vertices.back()(1)<<endl;
                }
            }

            auto c = cells_.emplace(
                Cell::id(real_i + real_i2, real_j + real_j2, layer_id),
                Cell(move(position), move(vertices), orientation, real_i + real_i2, real_j + real_j2, layer_id)
            );
            if(!c.second) {
                cout << "Warning: Cell with indices" << real_i + real_i2<< " " <<
                            real_j +real_j2 << ", "<<layer_id<<
                            " already exists (id="<<Cell::id(real_i + real_i2,real_j2+real_j, layer_id) << ")\n";
                cout << "Warning: This may indicate a bug in the id definition\n";
            }
            real_j2++;
        }
        real_i2++;
    }

    // // build histogram of cells for the selected layer
    if (display_layer == layer_id) {
        cell_histogram_.reset(new TH2Poly("cells", "cells",
            x_min>0?x_min*0.9:x_min*1.1, x_max>0?x_max*1.1:x_max*0.9,
            y_min>0?y_min*0.9:y_min*1.1, y_max>0?y_max*1.1:y_max*0.9));

        // take ownership of the histogram
        cell_histogram_->SetDirectory(0);

        for (auto& c : cells_) {

            auto& cell = c.second;
            if (cell.getLayer() == layer_id) {

                vector<double> binsx, binsy;
                for (const auto& vertex : cell.getVertices()) {
                    binsx.emplace_back(vertex(0));
                    binsy.emplace_back(vertex(1));
                }
                binsx.emplace_back(cell.getVertices()[0](0));
                binsy.emplace_back(cell.getVertices()[0](1));
                cell_histogram_->AddBin(binsx.size(), binsx.data(), binsy.data());
            }
        }
    }
}


void Geometry::print() {
    cout << "Printing the geometry : " << endl;
    if (klayer_ != -1)
        cout<<"the layer plane is "<<klayer_<< " at z position "<<parameters_.layers_z[klayer_]<<endl;
    else
        cout << "All the layer plane are simulated at the corresponding z position " << endl;
    for (auto& cell : cells_) {
        cout << "new cell with indices " << 
        "("<< cell.second.getIIndex() << ", " << cell.second.getJIndex() << ", " << cell.second.getLayer() << ")" <<
          " and position " << 
          "("<<cell.second.getX() << ", " << cell.second.getY() << ", "<< cell.second.getLayer() << ")" << endl;
    }
}

bool Geometry::isInCell(const TVectorD& position, const Cell& cell) const {
    // implementation below works for any convex cell described by its vertices
    // assumes vertices are consecutive along the cell perimeter and ordered along direct rotation

    // loop on pair of consective vertices 
    const auto& vertices = cell.getVertices();
    for (unsigned int i=0;i<vertices.size()-1;i++) {
        double xa = vertices[i](0);
        double ya = vertices[i](1);
        double xb = vertices[i+1](0);
        double yb = vertices[i+1](1);
        double sign = (xb-xa)*(position(1)-ya) - (yb-ya)*(position(0)-xa);
        if (sign<0.)
            return false;
    }
    return true;
}


// For the external geometry
const Cell* Geometry::closestCell(double x, double y) const {

  // beware that this function does not explicit check that the point is within
  // the cell, therefore if the cell grid is too small it will return
  // the closest cell (as the name indicates)

  // here we implement the a temporary bruteforce function for the more general geometry read from the json file
  // this needs to be optimize for timing efficiency

  // to test time spent here
  //return *cells_.begin();

    double r2min = std::numeric_limits<double>::max();
    auto icfound = cells_.cend();
    for (auto ic = cells_.cbegin(); ic != cells_.cend(); ic++) {
        double xcell = ic->second.getX();
        double ycell = ic->second.getY();
        // now compute the distances and take the smallest one
        double r2 = (xcell-x)*(xcell-x) + (ycell-y)*(ycell-y);
        if (r2 < r2min) {
            r2min = r2;
            icfound = ic;
        }
    }

    if (icfound == cells_.cend()) {
        // FIXME: do we want to throw an exception?
        throw std::string("Didn't find closest cell");
        //std::cout << "[Geometry::closestCell] Cell not found!! x, y " << 
          //x << " " << y << std::endl;
    }

  return &icfound->second;
}
