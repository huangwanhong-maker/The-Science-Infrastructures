// What a plug-in script receives as `api` from its create(api): the base
// classes, the result makers and the numerical helpers of ./plugin.js. A
// script written against this needs no import of the atlas; it may still
// import the domains and foundations by URL for the science it uses.
export * from './plugin.js';
