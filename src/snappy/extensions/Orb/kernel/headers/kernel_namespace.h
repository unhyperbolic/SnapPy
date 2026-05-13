/**
 *  @file kernel_namespace.h
 *  @brief Define namespace scope
 *
 *  If you wish to build the snappea kernel within a C++ namespace you
 *  can declare the namespace block in this file.  This is also a
 *  convenient place to put diagnostic pragmas.
 */
#ifdef _MSC_VER
#pragma warning(disable: 4190 4996)
#endif
#ifdef __APPLE__
#ifdef __cplusplus
#ifdef __clang__
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#endif
#endif
#endif

#ifndef SNAPPEA_NAMESPACE_BEGIN_SCOPE
// ifdef __cplusplus
// #define SNAPPEA_NAMESPACE_BEGIN_SCOPE namespace snappea {
// #else
#define SNAPPEA_NAMESPACE_BEGIN_SCOPE
// #endif
#endif

#ifndef SNAPPEA_NAMESPACE_END_SCOPE
// ifdef __cplusplus
// #define SNAPPEA_NAMESPACE_END_SCOPE }
// #else
#define SNAPPEA_NAMESPACE_END_SCOPE
// #endif
#define SNAPPEA_NAMESPACE_END_SCOPE
#endif
