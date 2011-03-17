/*!
 * \file
 * \author Chris de Pujo
 */

#include "MainPanel.hpp"

#include <algorithm>
#include <functional>

#include <wx/wx.h>

#include "../ConnectedComponents/trunk/cad.hpp"

BEGIN_EVENT_TABLE(MainPanel, wxPanel)
    EVT_PAINT(MainPanel::OnPaint)
END_EVENT_TABLE()

MainPanel::MainPanel(wxFrame * parent)
    : wxPanel(parent), oldWidth(-1), oldHeight(-1), nextLine(0),
      bMoreLinesToDraw(true)
{
    //enteredPolys.push_back("(x^2 + y^2)^4 - (x^2 - y^2)^2");
    //enteredPolys.push_back("x + y + 1/2");
    enteredPolys.push_back("x^2+y^2-1");
    enteredPolys.push_back("2*x+3*y^2-2");
    //enteredPolys.push_back("x^2+y^2-4");
    //enteredPolys.push_back("x+y");

    this->SetBackgroundColour(*wxWHITE);
}

void MainPanel::OnPaint(wxPaintEvent & event)
{
    wxPaintDC dc(this);

    this->nextLine = 0;
    this->bMoreLinesToDraw = true;
    render(dc);
}

void MainPanel::paintNow()
{
    wxClientDC dc(this);

    render(dc);
}
// bRestart is for OnPaint to force rendering to start over.
// Maybe just reset nextLine to zero?

// Also, maybe no need for bNewPolys? Just set old dimensions to -1?
void MainPanel::render(wxDC & dc)
{
    wxCoord width, height;
    dc.GetSize(&width, &height);

    this->OnPaintUpdatePolynomials(dc);

    dc.SetPen(*wxBLACK_PEN);

    static std::vector<int>     roots;
    const GiNaC::numeric        vBound(height);
    wxCoord                     vCoord = this->nextLine; // just a typedef'd int
    const long                  msDelta = 1000/8;

    assert(this->nextLine < height);

    wxStopWatch sw; // Maybe move this to the beginning of the function.

    for (GiNaC::numeric v = this->nextLine; v < vBound; ++v)
    {
        for (PolynomialQQ::vector::const_iterator
             i = this->cad->GetIrreducibles().begin(),
             e = this->cad->GetIrreducibles().end(); i != e; ++i)
            i->suby(v).addRoundedRootsTo(roots);

        std::sort(roots.begin(), roots.end());
        std::vector<int>::iterator i, e, newend =
            std::unique(roots.begin(), roots.end());
        roots.resize(newend - roots.begin());
        i = std::find_if(roots.begin(), roots.end(),
                         std::bind2nd(std::greater_equal<int>(), 0));
        e = std::find_if(roots.rbegin(), roots.rend(),
                         std::bind2nd(std::less<int>(), width)).base();

        for ( ; i != e; ++i)
        {
            /*if ((i+1) != e)
            {
                if (*(i+1) - *i > 1)
                {
                    dc.SetPen(this->colors[])
                    dc.DrawLine(*i + 1, vCoord, *(i+1), vCoord);
                }
            }*/
            dc.DrawPoint(*i, vCoord);
        }

        roots.clear();
        ++vCoord; // ghetto, but faster than converting.
        this->nextLine = vCoord;

        if (vCoord < height)
            if (sw.Time() >= msDelta)
                return;
    }

    this->bMoreLinesToDraw = false;
    this->nextLine = 0;
}

void MainPanel::OnPaintUpdatePolynomials(const wxDC & dc)
{
    wxCoord width, height;
    dc.GetSize(&width, &height);

    if (this->oldWidth == width && this->oldHeight == height)
        return;

    this->oldWidth = width;
    this->oldHeight = height;

    const GiNaC::numeric originx(width-1, 2);
    const GiNaC::numeric originy(height-1, 2);

    const GiNaC::numeric    ratio(3,4);
    GiNaC::numeric          factor;

    //

    PolynomialQ::vector P = CAD::Project(this->enteredPolys);
    Algebraic::vector points = PolynomialQ::FindRoots(P);
    P.clear();
    BOOST_FOREACH(Algebraic & alpha, points)
        alpha.takeAbs();
    Algebraic maxX = *std::max_element(points.begin(), points.end());
    points.clear();

    //PolynomialQQ::vector G(this->enteredPolys);
    PolynomialQQ::vector G(this->enteredPolys.begin(), this->enteredPolys.end());
    for (PolynomialQQ::vector::iterator i = G.begin(), e = G.end(); i != e; ++i)
        i->switchVariables();
    P = CAD::Project(G);
    points = PolynomialQ::FindRoots(P);
    P.clear();
    BOOST_FOREACH(Algebraic & alpha, points)
        alpha.takeAbs();
    Algebraic maxY = *std::max_element(points.begin(), points.end());
    points.clear();

    maxX.mul(GiNaC::numeric(width/2).inverse());      // maxX / width
    maxY.mul(GiNaC::numeric(height/2).inverse());     // maxY / height

    if (maxX >= maxY)
        factor = maxX.Approximate()/ratio;
    else
        factor = maxY.Approximate()/ratio;

    //

    /*PolynomialQQ::vector temp;
    this->irreducibles.clear();
    BOOST_FOREACH(const PolynomialQQ & f, this->F)
        f.addIrreducibleFactorsTo(temp);
    for (PolynomialQQ::vector::iterator i = temp.begin(); i != temp.end(); ++i)
        if (find(this->irreducibles.begin(), this->irreducibles.end(), *i) ==
            this->irreducibles.end())
            this->irreducibles.push_back(*i);
    BOOST_FOREACH(PolynomialQQ & f, this->irreducibles)
        f.linearSubs(factor,        originx*(-1)*factor,
                     factor*(-1),   originy*factor);*/

    PolynomialQQ::vector temp(this->enteredPolys);
    for (PolynomialQQ::vector::iterator i = temp.begin(), e = temp.end();
         i != e; ++i)
        i->linearSubs(factor,       originx*(-1)*factor,
                      factor*(-1),  originy*factor);

    cad = std::auto_ptr<CAD>(new CAD(temp));
}

/*void OnPress() {
    this->oldWidth = this->oldHeight = -1;
    bMoreLinesToDraw = true;
    paintNow();
}*/
