# Build xsd

list(APPEND projects xsd)
set(xsd_url "http://www.codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2")

ExternalProject_Add(xsd
    URL ${xsd_url}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(xsd source_dir)
set(XSD_DIR "${source_dir}" CACHE INTERNAL "")
