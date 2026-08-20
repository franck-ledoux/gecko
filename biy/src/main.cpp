// biy ("Block It Yourself"): loads a .msh file as a gecko GeomModel, shows it in an interactive 3D
// window, and lets you build a Blocking against it — from the on-screen panel, from mouse drags on
// the block's corners, or from the Python console on stdin. See docs/user-guide/biy.md.
#include <cstdlib>
#include <exception>
#include <iostream>

#include <polyscope/polyscope.h>

#include "BiyApp.h"
#include "PythonConsole.h"

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: biy <model.msh>\n";
        return 1;
    }

    try {
        polyscope::init();

        gecko::biy::BiyApp app(argv[1]);
        polyscope::state::userCallback = [&app] { app.per_frame(); };

        gecko::biy::PythonConsole console(app);

        // Verification aid: render a few frames, save what the window shows, and exit — lets the
        // rendering path be checked without a human at the keyboard.
        if (const char *shot = std::getenv("BIY_SCREENSHOT")) {
            for (int i = 0; i < 5; ++i) polyscope::frameTick();
            polyscope::screenshot(shot);
            console.stop();
            std::cout << "Wrote screenshot to " << shot << " — press Enter to exit.\n" << std::flush;
            return 0;
        }

        // frameTick() instead of show(): it renders one frame and returns, leaving this loop free to
        // also watch for the console asking to quit.
        while (!polyscope::windowRequestsClose() && console.running()) {
            polyscope::frameTick();
        }

        if (console.running()) {
            std::cout << "\nWindow closed — press Enter to exit.\n" << std::flush;
            console.stop();
        }
    } catch (const std::exception &e) {
        std::cerr << "biy: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
