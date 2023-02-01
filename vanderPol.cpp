#include <iostream>
#include <iomanip>
#include <string>
#include "dg/algorithm.h"
#include "dg/file/file.h"

namespace equations
{

double l2norm( const std::array<double,2>& y)
{
    return sqrt( y[0]*y[0]+y[1]*y[1]);
}

struct Implicit;
struct ImplicitSolver;

struct Explicit
{
    friend class Implicit;
    friend class ImplicitSolver;
    Explicit(){};
    Explicit( Json::Value js, bool partitioned)
    {
        m_equations = js["physical"].get( "equations", "van-der-Pol").asString();
        //This is the stiffness parameter
        m_sigma = js["physical"].get( "stiffness", 5.).asDouble();
        m_partitioned = partitioned;
    }
    void operator()( double t, const std::array<double,2>& y, std::array<double,2>& yp)
    {
        m_counter++;
        if( m_equations == "van-der-Pol")
        {
            yp[0] = y[1];
            yp[1] = - y[0];
            if( !m_partitioned)
                yp[1] += m_sigma*(1.-y[0]*y[0])*y[1];
        }
        else if ( m_equations == "Pareschi-Russo")
        {
            yp[0] = -y[1];
            yp[1] = y[0];
            if( !m_partitioned)
                yp[1] += m_sigma*(sin(y[0])-y[1]);
        }
    }
    unsigned get_counter() const{return m_counter;}
    private:
    std::string m_equations;
    double m_sigma;
    bool m_partitioned;
    unsigned m_counter = 0.;
};

struct Implicit
{
    Implicit(){}
    Implicit( Explicit& ex):m_ex(&ex){}
    void operator()( double t, const std::array<double,2>& y, std::array<double,2>& yp)
    {
        m_counter++;
        yp[0] = yp[1] = 0.;
        if( m_ex->m_equations == "van-der-Pol")
        {
            if( m_ex->m_partitioned)
                yp[1] = m_ex->m_sigma*(1.-y[0]*y[0])*y[1];
        }
        else if ( m_ex->m_equations == "Pareschi-Russo")
        {
            if( m_ex->m_partitioned)
                yp[1] = m_ex->m_sigma*(sin(y[0])-y[1]);
        }

    }
    unsigned get_counter() const{return m_counter;}
    private:
    Explicit* m_ex;
    unsigned m_counter = 0.;
};

static unsigned solve_counter = 0;
static unsigned iterations_counter = 0;
struct ImplicitSolver
{
    ImplicitSolver(){};
    ImplicitSolver( Explicit& ex, double eps_time):m_ex(&ex){
        m_eps_time = eps_time;
    }
    std::array<double,2> copyable(){ return {0,0};}
    // y + alpha I(t,y) = rhs
    void operator()( double alpha, double t, std::array<double,2>& y, const std::array<double,2>& rhs)
    {
        solve_counter++;
        y[0] = rhs[0];
        y[1] = rhs[1];

        if( m_ex->m_equations == "van-der-Pol")
        {
            if( m_ex->m_partitioned)
            {
                double prefactor = 1.+alpha*m_ex->m_sigma*(1.-y[0]*y[0]);
                if(prefactor == 0)
                    throw std::runtime_error( "Prefactor 0 in implicit solve");
                y[1] = rhs[1]/prefactor;
            }
            else
            {
                double eps = 1e100;
                iterations_counter = 0; // avoid deadlock
                //std::vector<double> convrate;
                while( eps > m_eps_time && iterations_counter < 100)
                {
                    double fx = y[0]+alpha*y[1]-rhs[0];
                    double fy = y[1]+alpha*m_ex->m_sigma*y[1]*(1.-y[0]*y[0])
                        - alpha*y[0]-rhs[1];
                    double Jxx = 1., Jxy = alpha;
                    double Jyx = -2.*alpha*m_ex->m_sigma*y[0]*y[1]-alpha;
                    double Jyy = 1.+alpha*m_ex->m_sigma*(1.-y[0]*y[0]);
                    double det = Jxx*Jyy-Jxy*Jyx;
                    if(det == 0)
                        throw std::runtime_error( "Determinant 0 in implicit solve");
                    double y_old0 = y[0], y_old1 = y[1];
                    y[0] = y[0] - 1./det*(Jyy*fx - Jxy*fy);
                    y[1] = y[1] - 1./det*(-Jyx*fx + Jxx*fy);
                    double eps_old = eps;
                    eps = sqrt( (y[0]-y_old0)*(y[0]-y_old0) + (y[1]-y_old1)*(y[1]-y_old1));
                    //convrate.append( eps/eps_old);
                    //std::cout << iterations_counter << " "<<alpha<< " " <<eps<<" "<<eps/eps_old<< "\n";
                    if( eps > 10*eps_old) return; // bad convergence
                    iterations_counter ++;
                }
            }
        }
        else if ( m_ex->m_equations == "Pareschi-Russo")
        {
            if( m_ex->m_partitioned)
            {
                double prefactor = 1.-alpha*m_ex->m_sigma;
                if(prefactor == 0)
                    throw std::runtime_error( "Prefactor 0 in implicit solve");
                y[1] = (rhs[1]-alpha*m_ex->m_sigma*sin(y[0]))/prefactor;
            }
            else
            {
                double eps = 1;
                iterations_counter = 0; // avoid deadlock
                while( eps > m_eps_time && iterations_counter < 100)
                {
                    double fx = y[0]-alpha*y[1]-rhs[0];
                    double fy = y[1]+alpha*m_ex->m_sigma*y[1]*(sin(y[0])-y[1])
                        + alpha*y[0]-rhs[1];
                    double Jxx = 1., Jxy = -alpha;
                    double Jyx = +alpha*m_ex->m_sigma*cos(y[0])+alpha;
                    double Jyy = 1.-alpha*m_ex->m_sigma;
                    double det = Jxx*Jyy-Jxy*Jyx;
                    if(det == 0)
                        throw std::runtime_error( "Determinant 0 in implicit solve");
                    double y_old0 = y[0], y_old1 = y[1];
                    y[0] = y[0] - 1./det*(Jyy*fx - Jxy*fy);
                    y[1] = y[1] - 1./det*(-Jyx*fx + Jxx*fy);
                    eps = sqrt( (y[0]-y_old0)*(y[0]-y_old0) + (y[1]-y_old1)*(y[1]-y_old1));
                    iterations_counter ++;
                }
            }
        }

    }
    private:
    Explicit* m_ex;
    double m_eps_time;
};


struct Variables{
    const std::array<double,2>& y0;
    const double& time;
    const double& dt;
    Json::Value& js;
    double eps0;
    unsigned nfailed, nsteps, nex, nim, nsolve, niter;
};

struct Record1d{
    std::string name;
    std::string long_name;
    std::function<double( Variables&)> function;
};

std::vector<Record1d> diagnostics1d_list = {
    {"first", "First variable y0",
        [](Variables& v){
            return v.y0[0];
        }
    },
    {"second", "Second variable y0",
        [](Variables& v){
            return v.y0[1];
        }
    },
    {"dt", "Proposed timestep",
        [](Variables& v){
            return v.dt;
        }
    },
    {"error", "Normalized error",
        [](Variables& v){
            return v.eps0;
        }
    },
    {"nfailed", "Number of failed steps",
        []( Variables& v ) {
            return v.nfailed;
        }
    },
    {"nsteps", "Number of calls to the timestepper (including failed steps)",
        [](Variables& v) {
            return v.nsteps;
        }
    },
    {"nex", "Number of calls to the explicit rhs",
        [](Variables& v) {
            return v.nex;
        }
    },
    {"nim", "Number of calls to the implicit rhs",
        [](Variables& v) {
            return v.nim;
        }
    },
    {"nsolve", "Number of calls to the implicit solver",
        [](Variables& v) {
            return v.nsolve;
        }
    },
    {"niter", "Number of iterations in the implicit solver",
        [](Variables& v) {
            return v.niter;
        }
    },
};
} //namespace equations

int main( int argc, char* argv[])
{
    ////Parameter initialisation ////////////////////////////////////////////
    Json::Value js;
    if( argc == 1)
        dg::file::file2Json( "input/default.json", js, dg::file::comments::are_discarded);
    else
        dg::file::file2Json( argv[1], js);
    std::cout << js <<std::endl;

    /////////////////////////////////////////////////////////////////
    std::array<double,2> y0 = {2.,0.};
    std::string equations = js["physical"].get( "equations", "van-der-Pol").asString();
    if( equations == "Pareschi-Russo")
        y0 = {M_PI/2., 1./2.};
    else if( !(equations == "van-der-Pol"))
        throw std::runtime_error( "Equations "+equations+" not recognized!\n");

    std::string type = js["timestepper"].get( "type", "adaptive").asString();
    double reject_limit = js["timestepper"].get( "reject-limit", 2).asDouble();
    // adaptive, adaptive-imex or adaptive-implicit
    std::string tableau = js["timestepper"].get("tableau",
            "Bogacki-Shampine-4-2-3").asString();
    double rtol= js["timestepper"].get("rtol", 1e-4).asDouble();
    double atol= js["timestepper"].get("atol", 1e-6).asDouble();
    double tend = js["output"].get( "tend", 1.0).asDouble();
    bool partitioned = false;
    if( type == "adaptive-imex")
        partitioned = "true";
    equations::Explicit ex( js, partitioned);
    equations::Implicit im( ex);
    double eps_time = 1.;
    if( type == "adaptive-implicit")
        eps_time = js["timestepper"].get( "eps_time", 1e-10).asDouble();
    equations::ImplicitSolver solver( ex, eps_time);

    dg::Adaptive<dg::ERKStep< std::array<double,2> >> adaptive;
    dg::Adaptive<dg::ARKStep< std::array<double,2>>> adaptive_imex;
    dg::Adaptive<dg::DIRKStep<std::array<double,2>>> adaptive_implicit;
    if( type == "adaptive")
        adaptive.construct( tableau, y0);
    else if( type == "adaptive-imex")
        adaptive_imex.construct( tableau, y0);
    else if( type == "adaptive-implicit")
        adaptive_implicit.construct( tableau, y0);
    else
        throw std::runtime_error( "Timestepper type "+type+" not recognized\n");
    std::string controller = js["timestepper"].get( "controller", "pid-control").asString();
    std::function< double( std::array<double,3>,std::array<double,3>,unsigned,unsigned)> control;
    if( controller == "i-control")
        control = dg::i_control;
    else if( controller == "pi-control")
        control = dg::pi_control;
    else if( controller == "pid-control")
        control = dg::pid_control;
    else if( controller == "ex-control")
        control = dg::ex_control;
    else if( controller == "im-control")
        control = dg::im_control;
    else if( controller == "imex-control")
        control = dg::imex_control;
    else
        throw std::runtime_error( "Controller "+controller+" not recognized\n");

    double dt = 1e-3, time = 0.;
    // Set up netcdf
    std::string inputfile = js.toStyledString(); //save input without comments, which is important if netcdf file is later read by another parser
    std::string outputfile;
    if( argc == 1 || argc == 2)
        outputfile = "vanDerPol.nc";
    else
        outputfile = argv[2];
    /// //////////////////////set up netcdf/////////////////////////////////////
    dg::file::NC_Error_Handle err;
    int ncid=-1;
    try{
        err = nc_create( outputfile.c_str(),NC_NETCDF4|NC_CLOBBER, &ncid);
    }catch( std::exception& e)
    {
        std::cerr << "ERROR creating file "<<outputfile<<std::endl;
        std::cerr << e.what()<<std::endl;
        return -1;
    }
    /// Set global attributes
    std::map<std::string, std::string> att;
    att["title"] = "Output file of timesteppers/vanDerPol.cpp";
    ///Get local time and begin file history
    auto ttt = std::time(nullptr);
    auto tm = *std::localtime(&ttt);

    std::ostringstream oss;
    ///time string  + program-name + args
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    for( int i=0; i<argc; i++) oss << " "<<argv[i];
    att["history"] = oss.str();
    att["comment"] = "Find more info in timesteppers/controller.ipynb";
    att["source"] = "FELTOR";
    att["references"] = "https://github.com/feltor-dev/feltor";
    att["inputfile"] = inputfile;
    for( auto pair : att)
        err = nc_put_att_text( ncid, NC_GLOBAL,
            pair.first.data(), pair.second.size(), pair.second.data());

    int dim_id, tvarID;
    std::map<std::string, int> id1d;
    err = dg::file::define_time( ncid, "time", &dim_id, &tvarID);
    //Create field IDs
    for( auto& record : equations::diagnostics1d_list)
    {
        std::string name = record.name;
        std::string long_name = record.long_name;
        id1d[name] = 0;
        err = nc_def_var( ncid, name.data(), NC_DOUBLE, 1, &dim_id,
                &id1d.at(name));
        err = nc_put_att_text( ncid, id1d.at(name), "long_name", long_name.size(),
            long_name.data());
    }
    err = nc_enddef(ncid);
    size_t start[1] = {0};
    size_t count[1] = {1};
    equations::Variables var = {y0, time, dt, js,1., 0,0,0,0,0,0};
    for( auto& record : equations::diagnostics1d_list)
    {
        double result = record.function( var);
        nc_put_vara_double( ncid, id1d.at(record.name), start, count, &result);
    }
    err = nc_put_vara_double( ncid, tvarID, start, count, &time);

    //////////////////Main time loop ///////////////////////////////////
    while( (tend - time) > 1e-14)
    {
        if( time+dt > tend)
            dt = tend-time;
        // Compute a step and error
        if( type == "adaptive")
        {
            adaptive.step( ex, time, y0, time, y0, dt,
                control, equations::l2norm, rtol, atol,reject_limit);
            var.eps0 = adaptive.get_error();
            if( adaptive.failed())
                var.nfailed++;
        }
        else if( type == "adaptive-imex")
        {
            adaptive_imex.step( std::tie(ex, im, solver), time, y0, time, y0, dt,
                control, equations::l2norm, rtol, atol,reject_limit);
            var.eps0 = adaptive_imex.get_error();
            if( adaptive_imex.failed())
                var.nfailed++;
        }
        else if( type == "adaptive-implicit")
        {
            adaptive_implicit.step( std::tie( im, solver), time, y0, time, y0, dt,
                control, equations::l2norm, rtol, atol,reject_limit);
            var.eps0 = adaptive_implicit.get_error();
            if( adaptive_implicit.failed())
                var.nfailed++;
        }
        var.nsteps++;
        var.nex = ex.get_counter();
        var.nim = im.get_counter();
        var.nsolve = equations::solve_counter;
        var.niter = equations::iterations_counter;
        /////////////////////////////////output/////////////////////////
        start[0]++;
        for( auto& record : equations::diagnostics1d_list)
        {
            double result = record.function( var);
            nc_put_vara_double( ncid, id1d.at(record.name), start, count, &result);
        }
        err = nc_put_vara_double( ncid, tvarID, start, count, &time);
    }
    err = nc_close(ncid);

    return 0;
}
