SET(DEAL_II_DIR "/home/liam/dev/deallii/installed" CACHE PATH "dealii directory")
FIND_PACKAGE(deal.II 8.2 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Include as SYSTEM to suppress "warning: extra ';' [-Wpedantic]" for deal.ii
INCLUDE_DIRECTORIES(SYSTEM ${DEAL_II_INCLUDE_DIRS})
