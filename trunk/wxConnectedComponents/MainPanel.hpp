/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __MAINPANEL__
#define __MAINPANEL__

#include <memory>

#include <wx/wx.h>

#include "../ConnectedComponents/trunk/PolynomialQQ.hpp"
#include "../ConnectedComponents/trunk/cad.hpp"

class MainPanel : public wxPanel
{
public:

    MainPanel(wxFrame * parent);

    void OnPaint(wxPaintEvent & event);
    void paintNow();
    void render(wxDC & dc);

    bool moreLinesToDraw() { return bMoreLinesToDraw; }

private:

    std::auto_ptr<CAD> cad;
    PolynomialQQ::vector enteredPolys;//, irreducibles;

    wxCoord oldWidth, oldHeight, nextLine;
    bool bMoreLinesToDraw;

    void OnPaintUpdatePolynomials(const wxDC & dc);

    DECLARE_EVENT_TABLE()
};

#endif // __MAINPANEL__
