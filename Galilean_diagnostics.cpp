// Diagnosticos locais: nao altera colisao, recoloracao, streaming ou biblioteca compartilhada.
// Populacoes pos-colisao: rho*u = soma_i f_i*c_i - F/2 (quando ini_force existe).
// Referencial principal: velocidade uniforme IMPOSTA v0 ex. Tambem mede no referencial Q/M.
#include <algorithm>
#include <limits>

struct GI_DATA {
    double mean_x=0, mean_y=0, mean_z=0, drop_x=0, drop_y=0;
    double max_ref=0, rms_ref=0, max_mean=0, rms_mean=0, max_lab=0, rms_lab=0;
    double min_f=std::numeric_limits<double>::infinity();
    double min_rho=std::numeric_limits<double>::infinity();
    double min_phi=1, max_phi=0, axis_ratio=1, center_dx=0, center_dy=0, x_unwrapped=0;
    long long negative_total=0, negative_color=0;
    int bad=0;
    std::vector<double> ux, uy, uz, phi;
};

inline double gi_image(double dx, double length) {
    return dx - length * std::floor(dx / length + 0.5);
}

// A medida de Laplace original ja usa centroide circular e distancias de imagem minima.
// Somente os indicadores secundarios de linha/canto precisam acompanhar o referencial periodico.
bool calc_galilean(GEOMETRY geometry, LATTICE lattice, DROP& drop, GI_DATA& g,
                   double v0, int step, double x_ini, double y_ini) {
    g = GI_DATA();
    const int n = geometry.fluid, nx=geometry.nx, ny=geometry.ny, nz=geometry.nz;
    if (nz != 1) { g.bad=1; return false; }
    g.ux.resize(n); g.uy.resize(n); g.uz.resize(n); g.phi.resize(n);
    double mass=0, mass_r=0, qx=0, qy=0, qz=0, qrx=0, qry=0;
    for (int p=0; p<n; ++p) {
        double rr=0, rb=0, mx=0, my=0, mz=0;
        for (int i=0; i<nvel; ++i) {
            const double fr=lattice.inif_R[p*nvel+i], fb=lattice.inif_B[p*nvel+i];
            const double f=fr+fb;
            if (!std::isfinite(fr) || !std::isfinite(fb)) ++g.bad;
            g.min_f=std::min(g.min_f,f);
            if (f < -1.e-12) ++g.negative_total;
            if (fr < -1.e-12) ++g.negative_color;
            if (fb < -1.e-12) ++g.negative_color;
            rr+=fr; rb+=fb;
            mx+=f*lattice.c_i[i*dim]; my+=f*lattice.c_i[i*dim+1]; mz+=f*lattice.c_i[i*dim+2];
        }
        const double rho=rr+rb;
        if (!std::isfinite(rho) || rho<=0 || !std::isfinite(mx) || !std::isfinite(my) || !std::isfinite(mz)) {
            ++g.bad; continue;
        }
        if (lattice.ini_force != nullptr) {
            const double* F=lattice.ini_force+3*p;
            if (!std::isfinite(F[0]) || !std::isfinite(F[1]) || !std::isfinite(F[2])) {
                ++g.bad; continue;
            }
            mx-=0.5*F[0]; my-=0.5*F[1]; mz-=0.5*F[2];
        }
        g.min_rho=std::min(g.min_rho,rho);
        g.ux[p]=mx/rho; g.uy[p]=my/rho; g.uz[p]=mz/rho; g.phi[p]=rr/rho;
        g.min_phi=std::min(g.min_phi,g.phi[p]); g.max_phi=std::max(g.max_phi,g.phi[p]);
        mass+=rho; mass_r+=rr; qx+=mx; qy+=my; qz+=mz;
        qrx+=rr*g.ux[p]; qry+=rr*g.uy[p];
    }
    if (g.bad || !(mass>0) || !(mass_r>0)) return false;
    g.mean_x=qx/mass; g.mean_y=qy/mass; g.mean_z=qz/mass;
    g.drop_x=qrx/mass_r; g.drop_y=qry/mass_r;
    calc_laplace_tension(geometry,lattice,drop);
    double sref=0, smean=0, slab=0, sxx=0, syy=0, sxy=0, sx=0, sy=0, weight=0;
    double left_w=0, right_w=0, left_x=0, right_x=0;
    const int y_line=((int)std::floor(drop.y0+0.5)) % ny;
    for (int pos=0; pos<nx*ny; ++pos) {
        if (!geometry.ini[pos]) continue;
        const int p=geometry.ini[pos]-1;
        const double ux=g.ux[p], uy=g.uy[p], uz=g.uz[p];
        const double a=(ux-v0)*(ux-v0)+uy*uy+uz*uz;
        const double b=(ux-g.mean_x)*(ux-g.mean_x)+(uy-g.mean_y)*(uy-g.mean_y)+(uz-g.mean_z)*(uz-g.mean_z);
        const double c=ux*ux+uy*uy+uz*uz;
        sref+=a; smean+=b; slab+=c;
        g.max_ref=std::max(g.max_ref,std::sqrt(a));
        g.max_mean=std::max(g.max_mean,std::sqrt(b));
        g.max_lab=std::max(g.max_lab,std::sqrt(c));
        const int x=pos%nx, y=pos/nx;
        const double dx=gi_image(x-drop.x0,nx), dy=gi_image(y-drop.y0,ny), phi=g.phi[p];
        weight+=phi; sx+=phi*dx; sy+=phi*dy;
        sxx+=phi*dx*dx; syy+=phi*dy*dy; sxy+=phi*dx*dy;
        if (y==y_line) {
            const double rr=density(lattice.inif_R+p*nvel), rb=density(lattice.inif_B+p*nvel), w=rr*rb;
            if (dx<0) {left_w+=w; left_x+=w*dx;} else {right_w+=w; right_x+=w*dx;}
        }
    }
    g.rms_ref=std::sqrt(sref/n); g.rms_mean=std::sqrt(smean/n); g.rms_lab=std::sqrt(slab/n);
    const double xx=sxx/weight-(sx/weight)*(sx/weight), yy=syy/weight-(sy/weight)*(sy/weight);
    const double xy=sxy/weight-(sx/weight)*(sy/weight), trace=xx+yy;
    const double disc=std::sqrt((xx-yy)*(xx-yy)+4*xy*xy);
    g.axis_ratio=(trace>disc) ? std::sqrt((trace+disc)/(trace-disc)) : std::numeric_limits<double>::infinity();
    g.center_dx=gi_image(drop.x0-(x_ini+v0*step),nx);
    g.center_dy=gi_image(drop.y0-y_ini,ny);
    g.x_unwrapped=x_ini+v0*step+g.center_dx;
    // Correcoes apenas dos diagnosticos secundarios; sigma principal e R_area*delta_p nao mudam.
    drop.raio_linha=(left_w>0 && right_w>0) ? 0.5*(right_x/right_w-left_x/left_w) : 0;
    const int corner_x=((int)std::floor(drop.x0+0.5)+nx/2)%nx;
    const int corner_y=((int)std::floor(drop.y0+0.5)+ny/2)%ny;
    const int pc=geometry.ini[corner_x+nx*corner_y]-1;
    drop.p_canto=density(lattice.inif_R+pc*nvel)*(drop.cs2_R>0 ? drop.cs2_R : lattice.c_s2)+density(lattice.inif_B+pc*nvel)*(drop.cs2_B>0 ? drop.cs2_B : lattice.c_s2);
    drop.sigma_linha=drop.raio_linha*(drop.p_centro-drop.p_canto);
    drop.u_max=g.max_ref; drop.u_rms=g.rms_ref;
    return std::isfinite(drop.sigma) && std::isfinite(g.max_ref) && std::isfinite(g.axis_ratio);
}

void gi_header(std::ostream& out) {
    out << "# step v0 Mach sigma R umax_ref urms_ref umax_mean urms_mean umax_lab mean_x mean_y mean_z "
           "drop_ux drop_uy x0 y0 x_unwrapped center_dx center_dy axis_ratio mass_R mass_B min_f "
           "negative_total negative_color min_rho min_phi max_phi qx qy qz p_in p_out delta_p urms_lab\n";
}

void gi_write(std::ostream& out, int step, double v0, double cs2, const DROP& d, const GI_DATA& g) {
    out << std::setprecision(17) << step << ' ' << v0 << ' ' << v0/std::sqrt(cs2) << ' '
        << d.sigma << ' ' << d.raio << ' ' << g.max_ref << ' ' << g.rms_ref << ' '
        << g.max_mean << ' ' << g.rms_mean << ' ' << g.max_lab << ' '
        << g.mean_x << ' ' << g.mean_y << ' ' << g.mean_z << ' ' << g.drop_x << ' ' << g.drop_y << ' '
        << d.x0 << ' ' << d.y0 << ' ' << g.x_unwrapped << ' ' << g.center_dx << ' ' << g.center_dy << ' '
        << g.axis_ratio << ' ' << d.massa_R << ' ' << d.massa_B << ' ' << g.min_f << ' '
        << g.negative_total << ' ' << g.negative_color << ' ' << g.min_rho << ' ' << g.min_phi << ' ' << g.max_phi << ' '
        << d.qx << ' ' << d.qy << ' ' << d.qz << ' ' << d.p_in << ' ' << d.p_out << ' ' << d.delta_p << ' ' << g.rms_lab << '\n';
    out.flush();
}

void gi_final_field(GEOMETRY geo, const GI_DATA& g, double v0) {
    std::ofstream out("galilean_final.vtk");
    out << std::setprecision(10) << "# vtk DataFile Version 3.0\nGalilean test; velocity relative to imposed v0\nASCII\n"
        << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << geo.nx << ' ' << geo.ny << " 1\nORIGIN 0 0 0\nSPACING 1 1 1\n"
        << "POINT_DATA " << geo.nx*geo.ny << "\nSCALARS phi double 1\nLOOKUP_TABLE default\n";
    for (int pos=0; pos<geo.nx*geo.ny; ++pos) out << g.phi[geo.ini[pos]-1] << '\n';
    out << "VECTORS velocity_relative double\n";
    for (int pos=0; pos<geo.nx*geo.ny; ++pos) {
        const int p=geo.ini[pos]-1;
        out << g.ux[p]-v0 << ' ' << g.uy[p] << ' ' << g.uz[p] << '\n';
    }
}
